/* ============================================================================
 * probfct_PMSM_speed_static.c
 * IPMSM SPEED-MPC problem definition, ported to the OLD (static/HLS) GRAMPC API
 * used by the crane2d synthesis framework.
 *
 * PORTED FROM: probfct_PMSM_speed.c (GRAMPC v2.3 API)
 * PORTED TO  : the older API in probfct.h of the static framework, i.e.
 *                - no  'const typeGRAMPCparam *param'  argument
 *                - xdes/udes passed explicitly to lfct/dldx/dldu
 *                - ffct has NO p argument
 *                - dfdx_vec/dfdu_vec arg order is (out,t,x,vec,u,p,userparam)
 *                  NOTE: 'vec' comes BEFORE 'u' (opposite of the new API!)
 *                - ocp_dim uses char* (not typeInt*)
 *
 * KEY DESIGN DECISION -- DISTURBANCE REMOVED (Np=0):
 *   The v2.3 file had Np=2 with p[0],p[1] as an additive voltage disturbance
 *   (offset-free). The old framework's ffct has no p argument, so the
 *   disturbance is DROPPED here (d = 0). This also makes the GRAMPC baseline a
 *   fair like-for-like comparison against the DNN clone, which likewise has no
 *   offset-free correction. If you later want the disturbance back, add d_d,d_q
 *   to userparam (see pSys map below) and add them in ffct's out[0],out[1].
 *
 * LUT HANDLING (unchanged in spirit from your v2.3 file):
 *   Ld and Lq are NOT looked up inside GRAMPC. They are frozen per control step
 *   by the external lut_ind block and passed in via userparam (pSys[1],pSys[2]).
 *   Likewise T_L comes from the external EKF via pSys[7].
 *
 * userparam LAYOUT (pSys):
 *   pSys[0] = R        stator resistance      (3.5)
 *   pSys[1] = Ld       <- FROM LUT, per step
 *   pSys[2] = Lq       <- FROM LUT, per step
 *   pSys[3] = psi_PM   flux linkage           (0.17)
 *   pSys[4] = zp       pole pairs             (3)
 *   pSys[5] = J        inertia                (4e-4)
 *   pSys[6] = mu       friction               (4e-4)
 *   pSys[7] = T_L      <- FROM EKF, per step
 *   pSys[8] = Umax^2   voltage circle bound   ((560/sqrt(3))^2)
 *   pSys[9] = Imax^2   current circle bound   (10^2)
 *   pSys[10..15] = pCost[0..5] = [q_id, q_iq, q_w, q_theta, r_ud, r_uq]
 *   (NPSYS = 10, so pCost = pSys + 10, exactly as the v2.3 file)
 *
 * STATE / CONTROL:
 *   x = [i_d, i_q, omega, theta]   (Nx=4)
 *   u = [u_d, u_q]                 (Nu=2)
 *   h = [voltage circle, current circle] <= 0   (Nh=2)
 * ==========================================================================*/

#include "probfct.h"

#define NPSYS 10

/* square macro */
#define POW2(a) ((a)*(a))

/** OCP dimensions -- OLD API uses char* **/
void ocp_dim(char *Nx, char *Nu, char *Np, char *Ng, char *Nh, char *NgT, char *NhT, typeUSERPARAM *userparam)
{
	*Nx  = 4;
	*Nu  = 2;
	*Np  = 0;    /* disturbance dropped: old ffct has no p argument (see header) */
	*Nh  = 2;    /* inequalities: voltage circle, current circle */
	*Ng  = 0;
	*NgT = 0;
	*NhT = 0;
}

/** System function f(t,x,u,userparam)  -- NOTE: no p in the old API
    ------------------------------------ **/
void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, typeUSERPARAM *userparam)
{
	ctypeRNum *pSys = (ctypeRNum*)userparam;
	/* pSys: 0=R 1=Ld 2=Lq 3=psi 4=zp 5=J 6=mu 7=TL */

	/* d-axis:  did/dt = (ud - R*id + Lq*w*iq)/Ld        (p[0] dropped) */
	out[0] = (u[0] - pSys[0] * x[0] + pSys[2] * x[1] * x[2]) / pSys[1];

	/* q-axis:  diq/dt = (uq - R*iq - (psi + Ld*id)*w)/Lq (p[1] dropped) */
	out[1] = (u[1] - pSys[0] * x[1] - (pSys[3] + pSys[1] * x[0])*x[2]) / pSys[2];

	/* mech:    dw/dt = 0.5*(-2*zp*TL + 3*zp^2*(psi + (Ld-Lq)*id)*iq - 2*mu*w)/J */
	out[2] = ((typeRNum)0.5*(-2 * pSys[4] * pSys[7] + 3 * POW2(pSys[4])*
		(pSys[3] + (pSys[1] - pSys[2])*x[0])*x[1] - 2 * pSys[6] * x[2])) / pSys[5];

	/* angle:   dtheta/dt = w */
	out[3] = x[2];
}

/** Jacobian df/dx multiplied by vec -- OLD arg order: (out,t,x,vec,u,p,userparam)
    !! 'vec' is the 4th arg here, 'u' the 5th (swapped vs the v2.3 API) !! **/
void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	ctypeRNum *pSys = (ctypeRNum*)userparam;

	out[0] = -((vec[0] * pSys[0]) / pSys[1]) + ((typeRNum)1.5*vec[2] * POW2(pSys[4])*(pSys[1] - pSys[2])*x[1]) /
		pSys[5] - (vec[1] * pSys[1] * x[2]) / pSys[2];
	out[1] = -((vec[1] * pSys[0]) / pSys[2]) + ((typeRNum)1.5*vec[2] * POW2(pSys[4])*(pSys[3] + (pSys[1] - pSys[2])*x[0])) /
		pSys[5] + (vec[0] * pSys[2] * x[2]) / pSys[1];
	out[2] = vec[3] - (vec[2] * pSys[6]) / pSys[5] - (vec[1] * (pSys[3] + pSys[1] * x[0])) /
		pSys[2] + (vec[0] * pSys[2] * x[1]) / pSys[1];
	out[3] = 0;
}

/** Jacobian df/du multiplied by vec -- OLD arg order: (out,t,x,vec,u,p,userparam) **/
void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	ctypeRNum *pSys = (ctypeRNum*)userparam;

	out[0] = vec[0] / pSys[1];    /* d/d ud : enters id eq /Ld */
	out[1] = vec[1] / pSys[2];    /* d/d uq : enters iq eq /Lq */
}

/** Jacobian df/dp multiplied by vec -- Np=0, so empty **/
void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}


/** Integral cost l(t,x,u,p,xdes,udes,userparam)  -- xdes/udes explicit in old API
    -------------------------------------------------- **/
void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum *pSys  = (ctypeRNum*)userparam;
	ctypeRNum *pCost = pSys + NPSYS;

	/* SPEED-MPC cost (unified speed regulation, MTPA-emergent):
	   l = q_i*(id^2 + iq^2) + q_w*(w - w*)^2 + r*(u - u*)^2
	   Current terms penalize toward ZERO so the minimum-current split that
	   delivers the speed-tracking torque emerges (MTPA from loss minimization).
	   xdes[2] = speed reference w*.  T_L enters the dynamics via pSys[7]. */
	out[0] = pCost[0] * POW2(x[0])
		+ pCost[1] * POW2(x[1])
		+ pCost[2] * POW2(x[2] - xdes[2])
		+ pCost[3] * POW2(x[3] - xdes[3])
		+ pCost[4] * POW2(u[0] - udes[0])
		+ pCost[5] * POW2(u[1] - udes[1]);
}

/** Gradient dl/dx **/
void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum *pSys  = (ctypeRNum*)userparam;
	ctypeRNum *pCost = pSys + NPSYS;

	out[0] = 2 * pCost[0] * x[0];             /* d/d id : q_i*id^2       */
	out[1] = 2 * pCost[1] * x[1];             /* d/d iq : q_i*iq^2       */
	out[2] = 2 * pCost[2] * (x[2] - xdes[2]); /* d/d w  : q_w*(w-w*)^2   */
	out[3] = 2 * pCost[3] * (x[3] - xdes[3]); /* d/d th : q_theta (=0)   */
}

/** Gradient dl/du **/
void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum *pSys  = (ctypeRNum*)userparam;
	ctypeRNum *pCost = pSys + NPSYS;

	out[0] = 2 * pCost[4] * (u[0] - udes[0]);
	out[1] = 2 * pCost[5] * (u[1] - udes[1]);
}

/** Gradient dl/dp -- Np=0 **/
void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
}


/** Terminal cost V(T,x(T),p,xdes,userparam) -- not used (TerminalCost off)
    ---------------------------------------- **/
void Vfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
}
void dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
}
void dVdp(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
}
void dVdT(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
}


/** Equality constraints g(...) = 0  -- Ng=0, none
    --------------------------------------------------- **/
void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
void dgdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Inequality constraints h(t,x,u,p,userparam) <= 0
    h[0] : voltage circle   (ud^2 + uq^2 - Umax^2)/Umax^2 <= 0
    h[1] : current circle   (id^2 + iq^2 - Imax^2)/Imax^2 <= 0
    (normalized by the bound, as in the v2.3 file -- keeps the multipliers scaled)
    ------------------------------------------------------ **/
void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	ctypeRNum *pSys = (ctypeRNum*)userparam;

	out[0] = ((POW2(u[0]) + POW2(u[1])) - pSys[8]) / pSys[8];   /* voltage circle */
	out[1] = ((POW2(x[0]) + POW2(x[1])) - pSys[9]) / pSys[9];   /* current circle */
}

/** Jacobian dh/dx multiplied by vec -- OLD arg order: (out,t,x,u,p,vec,userparam) **/
void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
	ctypeRNum *pSys = (ctypeRNum*)userparam;

	out[0] = (2 * vec[1] * x[0]) / pSys[9];   /* only the current circle sees x */
	out[1] = (2 * vec[1] * x[1]) / pSys[9];
	out[2] = 0;
	out[3] = 0;
}

/** Jacobian dh/du multiplied by vec **/
void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
	ctypeRNum *pSys = (ctypeRNum*)userparam;

	out[0] = (2 * vec[0] * u[0]) / pSys[8];   /* only the voltage circle sees u */
	out[1] = (2 * vec[0] * u[1]) / pSys[8];
}

/** Jacobian dh/dp multiplied by vec -- Np=0 **/
void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Terminal equality constraints -- NgT=0, none
    -------------------------------------------------------- **/
void gTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
void dgTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
void dgTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
void dgTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Terminal inequality constraints -- NhT=0, none
    ----------------------------------------------------------- **/
void hTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
void dhTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
void dhTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
void dhTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Semi-implicit / RODAS functions -- not used (Integrator = Euler/Heun)
    ------------------------------------------------------- **/
void dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
void dfdxtrans(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
void dfdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
void dHdxdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *adj, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
void Mfct(typeRNum *out, typeUSERPARAM *userparam)
{
}
void Mtrans(typeRNum *out, typeUSERPARAM *userparam)
{
}
