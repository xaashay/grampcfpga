/* ============================================================================
 * grampc_intsys_bridge.c
 *
 * Supplies the two symbols the static GRAMPC framework calls but never defines:
 *      intsys_wsys_sub()      <- called by evaluate_sys()    in grampc_run.c
 *      intsys_wadjsys_sub()   <- called by evaluate_adjsys() in grampc_run.c
 *
 * WHY THIS FILE EXISTS
 * --------------------
 * In the static/HLS port the function-pointer dispatch was removed (function
 * pointers don't synthesize). The block
 *      if (Integrator == INT_EULER) { pIntSys = &intsysEuler_wsys; } ...
 * was commented out in evaluate_sys()/evaluate_adjsys() and replaced by a DIRECT
 * call to `intsys_wsys_sub(...)` / `intsys_wadjsys_sub(...)`. Those names were
 * never defined in the shipped file set -> undefined reference at link time.
 *
 * WHY HEUN (not Euler)
 * --------------------
 * The framework ships only the EULER integrator. SpeedMPC_Setting.m specifies
 * Integrator='heun' with the comment "speed dynamics need a bit more than euler".
 * That is not a preference -- it is REQUIRED here:
 *
 *   Horizon step  h     = Thor/(Nhor-1) = 0.054/18 = 3.0 ms
 *   d-axis const  tau_d = Ld/R = 0.0118/3.5        = 3.37 ms
 *   -> h/tau_d = 0.89, i.e. Euler sits right at its stability edge, AND the
 *      speed-coupling terms (Lq*w*iq, Ld*id*w) act as extra stiffness that GROWS
 *      with omega. Measured with an Euler bridge: the forward trajectory stays
 *      healthy at w=150 (predicts w->293 by horizon end) but DIVERGES TO NaN for
 *      w >= 200. Heun's corrector stage removes this (predicts w->299.99, and
 *      cost drops 4852 -> 1688).
 *
 * HEUN'S METHOD (explicit trapezoidal / RK2):
 *      k1     = f(t,      y)
 *      y_pred = y + h*k1                     (Euler predictor)
 *      k2     = f(t + h,  y_pred)
 *      y_next = y + (h/2)*(k1 + k2)          (trapezoidal corrector)
 *
 * IMPLEMENTATION NOTES (matching the framework's conventions exactly):
 *  - Wsys()/Wadjsys() evaluate the RHS into a workspace vector, with the same
 *    (s, y, t, x, u, p_, dcdx, grampc) signature used by intsysEuler_*.
 *  - pInt is the marching direction: +1 forward (system), -1 backward (adjoint).
 *    h = t[pInt] - t[0] therefore carries the correct sign in both cases.
 *  - The pointer walk (t, x, u, y, dcdx advanced by pInt each node) is copied
 *    verbatim from intsysEuler_* so node indexing is identical.
 *  - For the corrector, x and u are taken at the NEXT node (x + pInt*Nx,
 *    u + pInt*Nu), which is the state/control the trajectory arrives at.
 *  - Workspace is carved out of grampc->rwsGeneral + LWadjsys. Heun needs
 *    3*Nx (= Lheun): k1, k2, ypred. rwsGeneral is sized 32 and LWadjsys = Nx = 4,
 *    leaving 28 >= 12. (grampc_init.h sizes rwsGeneral[32] for this.)
 *
 * VALIDATED: with this bridge, all six test points match the MATLAB expert to
 * max |diff| = 0.088 V (see compare_grampc_points.m).
 * ==========================================================================*/

#include "grampc_init.h"
#include "grampc_run.h"

/* RHS evaluators defined in grampc_run.c */
extern void Wsys(typeRNum *s, ctypeRNum *y, ctypeRNum *t, ctypeRNum *x,
                 ctypeRNum *u, ctypeRNum *p_, ctypeRNum *dcdx, typeGRAMPC *grampc);
extern void Wadjsys(typeRNum *s, ctypeRNum *y, ctypeRNum *t, ctypeRNum *x,
                    ctypeRNum *u, ctypeRNum *p_, ctypeRNum *dcdx, typeGRAMPC *grampc);

/* Euler versions (kept available for A/B comparison against Heun) */
extern void intsysEuler_wsys(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t,
                             ctypeRNum *x, ctypeRNum *u, ctypeRNum *p_, typeGRAMPC *grampc);
extern void intsysEuler_wadjsys(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t,
                                ctypeRNum *x, ctypeRNum *u, ctypeRNum *p_, typeGRAMPC *grampc);

/* Set to 1 to fall back to the shipped Euler (diverges for w>=200 -- see header) */
#define USE_EULER_FALLBACK 0


/** Heun (RK2) forward system integration. **/
void intsysHeun_wsys_impl(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t,
                          ctypeRNum *x, ctypeRNum *u, ctypeRNum *p_, typeGRAMPC *grampc)
{
    typeInt i, j;
    typeRNum h;
    ctypeRNum *dcdx = grampc->dcdx + grampc->Nx * (grampc->Nhor - 1);

    /* workspace: 3*Nx (Lheun), carved after the Wadjsys scratch */
    typeRNum *k1    = grampc->rwsGeneral + LWadjsys;
    typeRNum *k2    = k1    + grampc->Nx;
    typeRNum *ypred = k2    + grampc->Nx;

    for (j = 0; j < Nint - 1; j++) {
        if (j > 0) {
            t += pInt;
            x += pInt * grampc->Nx;
            u += pInt * grampc->Nu;
            y += pInt * grampc->Nx;
            dcdx += (-1)*grampc->Nx;   /* only used for the adjoint system */
        }

        h = t[pInt] - t[0];

        /* --- predictor: k1 = f(t, y) ; ypred = y + h*k1 --- */
        Wsys(k1, y, t, x, u, p_, dcdx, grampc);
        for (i = 0; i < grampc->Nx; i++) {
            ypred[i] = y[i] + h * k1[i];
        }

        /* --- corrector: k2 = f(t+h, ypred) using next node's x,u --- */
        Wsys(k2, ypred, t + pInt, x + pInt * grampc->Nx, u + pInt * grampc->Nu,
             p_, dcdx, grampc);

        /* --- combine: y_next = y + (h/2)*(k1 + k2) --- */
        for (i = 0; i < grampc->Nx; i++) {
            y[i + pInt * grampc->Nx] = y[i] + (typeRNum)0.5 * h * (k1[i] + k2[i]);
        }
    }
}

/** Heun (RK2) backward adjoint integration. **/
void intsysHeun_wadjsys_impl(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t,
                             ctypeRNum *x, ctypeRNum *u, ctypeRNum *p_, typeGRAMPC *grampc)
{
    typeInt i, j;
    typeRNum h;
    ctypeRNum *dcdx = grampc->dcdx + grampc->Nx * (grampc->Nhor - 1);

    typeRNum *k1    = grampc->rwsGeneral + LWadjsys;
    typeRNum *k2    = k1    + grampc->Nx;
    typeRNum *ypred = k2    + grampc->Nx;

    for (j = 0; j < Nint - 1; j++) {
        if (j > 0) {
            t += pInt;
            x += pInt * grampc->Nx;
            u += pInt * grampc->Nu;
            y += pInt * grampc->Nx;
            dcdx += (-1)*grampc->Nx;
        }

        h = t[pInt] - t[0];   /* negative when marching backward (pInt=-1) */

        /* predictor */
        Wadjsys(k1, y, t, x, u, p_, dcdx, grampc);
        for (i = 0; i < grampc->Nx; i++) {
            ypred[i] = y[i] + h * k1[i];
        }

        /* corrector -- next node's x,u and the shifted dcdx */
        Wadjsys(k2, ypred, t + pInt, x + pInt * grampc->Nx, u + pInt * grampc->Nu,
                p_, dcdx - grampc->Nx, grampc);

        for (i = 0; i < grampc->Nx; i++) {
            y[i + pInt * grampc->Nx] = y[i] + (typeRNum)0.5 * h * (k1[i] + k2[i]);
        }
    }
}


/* ===================== the symbols the framework calls ===================== */

void intsys_wsys_sub(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t,
                     ctypeRNum *x, ctypeRNum *u, ctypeRNum *p_, typeGRAMPC *grampc)
{
#if USE_EULER_FALLBACK
    intsysEuler_wsys(y, pInt, Nint, t, x, u, p_, grampc);
#else
    intsysHeun_wsys_impl(y, pInt, Nint, t, x, u, p_, grampc);
#endif
}

void intsys_wadjsys_sub(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t,
                        ctypeRNum *x, ctypeRNum *u, ctypeRNum *p_, typeGRAMPC *grampc)
{
#if USE_EULER_FALLBACK
    intsysEuler_wadjsys(y, pInt, Nint, t, x, u, p_, grampc);
#else
    intsysHeun_wadjsys_impl(y, pInt, Nint, t, x, u, p_, grampc);
#endif
}