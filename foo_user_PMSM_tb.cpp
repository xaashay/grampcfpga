/* ============================================================================
 * foo_user_PMSM_tb.cpp -- C-simulation testbench for the GRAMPC IPMSM speed-MPC.
 *
 * Drives foo_user_PMSM at operating points spanning the speed range and checks
 * the returned [ud,uq] against the values verified in software, which were
 * themselves validated against the MATLAB expert (max |diff| = 0.088 V over all
 * six points -- see compare_grampc_points.m).
 *
 * WHY EACH SOLVE IS RESET
 * -----------------------
 * The static GRAMPC object is a GLOBAL, initialised once at load. After a
 * grampc_run its auglag multipliers / penalties / line-search state persist.
 * foo_user_PMSM_reset() restores them so every point is an independent COLD
 * solve -- which is what the reference values were produced with, and what
 * Controller_SpeedMPC does (it sets u0=[0,0] each step).
 *
 * TOLERANCE
 * ---------
 * Reference values are the double-precision software results. C-sim runs the
 * same double code, so agreement should be near-exact; TOL is loose enough to
 * absorb compiler/FP-reassociation differences but tight enough to catch a real
 * regression. If you later convert foo_data.h to ap_fixed, EXPECT this to fail
 * and re-baseline deliberately -- GRAMPC's gradient/line-search maths is far
 * more precision-sensitive than the DNN's feedforward pass.
 * ==========================================================================*/
#include <cstdio>
#include <cmath>
#include "foo_data.h"

extern void foo_user_PMSM(data_t_x_hat_in x_hat_in_int[X_HAT_IN_LENGTH],
                          data_t_u_opt_out u_opt_out_int[U_OPT_OUT_LENGTH]);
extern "C" void foo_user_PMSM_reset(void);

#define TOL      0.5      /* volts */
#define U_MAX    323.3161582

int main()
{
    /* [ id,  iq,   w,  theta,  Ld,     Lq,     T_L,  w_ref ] */
    static const double pts[6][8] = {
        {  0.0, 0.0,   0.0, 0.0, 0.0120, 0.0240, 3.0, 300.0 },
        { -0.5, 2.0,  50.0, 0.1, 0.0119, 0.0238, 3.0, 300.0 },
        { -1.0, 4.0, 150.0, 0.5, 0.0118, 0.0236, 3.0, 300.0 },
        { -1.2, 4.5, 290.0, 1.0, 0.0117, 0.0235, 3.0, 300.0 },
        { -1.2, 4.5, 300.0, 1.0, 0.0117, 0.0235, 3.0, 300.0 },
        { -1.2, 4.5, 300.0, 1.0, 0.0117, 0.0235, 6.0, 300.0 },
    };
    /* golden [ud,uq] -- software reference, validated vs MATLAB to 0.088 V */
    static const double gold[6][2] = {
        {  18.070, 323.316 },
        { -40.882, 163.484 },
        { -61.925,  95.735 },
        { -55.127,  48.699 },
        { -53.266,  44.598 },
        { -88.526,  69.595 },
    };
    static const char *lbl[6] = {
        "standstill->300", "accelerating", "mid-speed",
        "near target", "at target TL=3", "at target TL=6"
    };

    int fails = 0, nans = 0;
    double worst = 0.0;

    printf("=== GRAMPC IPMSM speed-MPC : C-simulation ===\n");
    printf("%-16s %-9s %-22s %-22s %8s\n", "point", "w [rad/s]", "got [ud,uq]", "golden [ud,uq]", "maxdiff");
    printf("--------------------------------------------------------------------------------------\n");

    for (int k = 0; k < 6; k++) {
        foo_user_PMSM_reset();                    /* independent cold solve */

        data_t_x_hat_in  in[X_HAT_IN_LENGTH];
        data_t_u_opt_out out[U_OPT_OUT_LENGTH];
        for (int i = 0; i < 8; i++) in[i] = (data_t_x_hat_in)pts[k][i];

        foo_user_PMSM(in, out);

        double ud = (double)out[0], uq = (double)out[1];
        bool bad = false;

        if (std::isnan(ud) || std::isnan(uq) || std::isinf(ud) || std::isinf(uq)) {
            nans++; bad = true;
        }
        double d0 = fabs(ud - gold[k][0]);
        double d1 = fabs(uq - gold[k][1]);
        double d  = (d0 > d1) ? d0 : d1;
        if (!bad && d > worst) worst = d;
        if (!bad && d > TOL) { fails++; bad = true; }

        printf("%-16s %9.1f [%8.3f,%8.3f]   [%8.3f,%8.3f]   %8.3f %s\n",
               lbl[k], pts[k][2], ud, uq, gold[k][0], gold[k][1], d,
               bad ? "<-- FAIL" : "");
    }

    /* --- envelope check: the box constraint umax must hold on each axis --- */
    printf("\n--- constraint check ---\n");
    int box_viol = 0;
    for (int k = 0; k < 6; k++) {
        foo_user_PMSM_reset();
        data_t_x_hat_in  in[X_HAT_IN_LENGTH];
        data_t_u_opt_out out[U_OPT_OUT_LENGTH];
        for (int i = 0; i < 8; i++) in[i] = (data_t_x_hat_in)pts[k][i];
        foo_user_PMSM(in, out);
        double ud = fabs((double)out[0]), uq = fabs((double)out[1]);
        if (ud > U_MAX + 1e-3 || uq > U_MAX + 1e-3) {
            box_viol++;
            printf("  point %d: box |u| exceeded (|ud|=%.3f |uq|=%.3f > %.3f)\n", k, ud, uq, U_MAX);
        }
    }
    if (!box_viol) printf("  box constraint |ud|,|uq| <= Umax: OK on all points\n");
    printf("  NOTE: the CIRCULAR limit (ud^2+uq^2 <= Umax^2) is a soft auglag\n");
    printf("        constraint -- small violations are normal for this real-time\n");
    printf("        suboptimal solver and are NOT a failure.\n");

    printf("\n=== RESULT ===\n");
    printf("  NaN/Inf   : %d\n", nans);
    printf("  mismatches: %d  (tol %.2f V)\n", fails, (double)TOL);
    printf("  worst diff: %.4f V\n", worst);
    printf("  %s\n", (nans == 0 && fails == 0) ? "PASS" : "FAIL");

    return (nans == 0 && fails == 0) ? 0 : 1;
}
