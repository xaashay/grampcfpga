/* ============================================================================
 * foo_user_PMSM.cpp -- IPMSM SPEED-MPC top for the static/HLS GRAMPC framework.
 *
 * ARCHITECTURE (identical to the DNN path -- LUT and EKF stay OUTSIDE GRAMPC):
 *     lut_ind(id,iq) -> Ld,Lq  ----\
 *     ekf_top(...)   -> T_L    -----+--> userparam --> grampc_run --> [ud,uq]
 *     w_ref                    ----/
 *   Ld,Lq,T_L are FROZEN per control step and written into userparam, exactly as
 *   Controller_SpeedMPC.m does (grampc.userparam(2)=Ld,(3)=Lq,(8)=TLhat).
 *
 * The static GRAMPC object itself lives in initialize_PMSM.h (framework pattern).
 * ==========================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include "foo_data.h"
#include "grampc_run.h"
#include "initialize_PMSM.h"

#define NX  4
#define NU  2

/* Reset the solver workspace to its initial state. The static object is a
 * global initialised once at load; after a grampc_run its auglag multipliers,
 * penalties and line-search state carry over. GRAMPC is designed to be called
 * repeatedly (warm-started) -- but the FIRST call works and later ones NaN,
 * which means some workspace is not being re-established. Call this before each
 * independent solve to reproduce a cold start. */
extern "C" void foo_user_PMSM_reset(void)
{
    typeInt i;
    for (i = 0; i < N_c * horizon; i++) {
        grampc->rws_pen[i]  = grampc->PenaltyMin;
        grampc->mult[i]     = 0;
        grampc->rws_cfct[i] = 0;
        grampc->cfctprev[i] = 0;
    }
    for (i = 0; i < N_x * horizon; i++) { grampc->x[i] = 0; grampc->adj[i] = 0; }
    for (i = 0; i < N_u * horizon; i++) {
        grampc->u[i] = 0; grampc->uls[i] = 0; grampc->uprev[i] = 0;
        grampc->gradu[i] = 0; grampc->graduprev[i] = 0; grampc->dcdu[i] = 0;
    }
    for (i = 0; i < N_x * (horizon + 1); i++) grampc->dcdx[i] = 0;
    grampc->lsExplicit[0] = (typeRNum)85.17561648;
    grampc->lsExplicit[1] = (typeRNum)157858.3299;
    grampc->lsExplicit[2] = (typeRNum)0.0001;
    grampc->lsExplicit[3] = (typeRNum)1;
}

extern "C" void foo_diag(void)
{
    printf("    x traj[0..3] : ");
    for(int i=0;i<4;i++) printf("%10.4f ", (double)grampc->x[i]);
    printf("\n    x traj[last] : ");
    int L=(horizon-1)*N_x;
    for(int i=0;i<4;i++) printf("%10.4f ", (double)grampc->x[L+i]);
    printf("\n    u traj[0..1] : %10.4f %10.4f\n",(double)grampc->u[0],(double)grampc->u[1]);
    printf("    J=[%g,%g] status=%d\n",(double)grampc->J[0],(double)grampc->J[1],(int)grampc->status);
}

void foo_user_PMSM( data_t_x_hat_in x_hat_in_int[X_HAT_IN_LENGTH],
                    data_t_u_opt_out u_opt_out_int[U_OPT_OUT_LENGTH])
{
    /* ---------- per-step inputs ---------- */
    /* state x0 = [id, iq, omega, theta] */
    for (int i = 0; i < NX; i++)
        grampc->x0[i] = (typeRNum)x_hat_in_int[i];

    /* Per-step params into userparam. NOTE: typeUSERPARAM is 'const typeRNum'
     * (the crane never wrote userparam -- its params were compile-time consts).
     * We MUST write Ld,Lq,T_L every step, and the storage IS real (an array
     * inside the struct), so cast away the const qualifier on write. */
    typeRNum *up = (typeRNum*)grampc->userparam;

    /* frozen Ld,Lq from the external LUT  -> pSys[1], pSys[2] */
    up[1] = (typeRNum)x_hat_in_int[IDX_LD];
    up[2] = (typeRNum)x_hat_in_int[IDX_LQ];

    /* load-torque estimate from the external EKF -> pSys[7] */
    up[7] = (typeRNum)x_hat_in_int[IDX_TL];

    /* speed reference -> xdes[2] */
    grampc->xdes[2] = (typeRNum)x_hat_in_int[IDX_WREF];

    /* cold start each solve (matches Controller_SpeedMPC: u0=[0,0]) */
    grampc->u0[0] = 0;  grampc->u0[1] = 0;

    /* ---------- solve ---------- */
    grampc_run(grampc);

    /* ---------- outputs ---------- */
    u_opt_out_int[0] = (data_t_u_opt_out)grampc->unext[0];
    u_opt_out_int[1] = (data_t_u_opt_out)grampc->unext[1];
}
