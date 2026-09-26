/* ============================================================================
 * foo_data.h -- I/O types and vector lengths for the IPMSM speed-MPC GRAMPC top.
 *
 * The crane framework's foo_user.cpp referenced this header but it wasn't in the
 * provided file set, so it's reconstructed here for the PMSM port.
 *
 * INPUT VECTOR (x_hat_in), length 8:
 *   [0] id      measured/estimated d-current      [A]
 *   [1] iq      measured/estimated q-current      [A]
 *   [2] omega   speed                             [rad/s elec]
 *   [3] theta   angle                             [rad elec]
 *   [4] Ld      frozen inductance   <- lut_ind    [H]
 *   [5] Lq      frozen inductance   <- lut_ind    [H]
 *   [6] T_L     load-torque estimate <- ekf_top   [Nm]
 *   [7] w_ref   speed reference (omega*)          [rad/s elec]
 *
 * OUTPUT VECTOR (u_opt_out), length 2:
 *   [0] u_d     d-axis voltage                    [V]
 *   [1] u_q     q-axis voltage                    [V]
 *
 * TYPES: double for the first (software / C-sim) pass, so results can be checked
 * bit-for-bit against MATLAB before any fixed-point conversion. For an HLS
 * fixed-point build, switch these to ap_fixed and include ap_fixed.h -- but do
 * that ONLY after the double version is verified against MATLAB, since GRAMPC's
 * gradient/line-search maths is far more sensitive to precision than the DNN's
 * feedforward pass.
 * ==========================================================================*/
#ifndef FOO_DATA_H
#define FOO_DATA_H

#define X_HAT_IN_LENGTH    8
#define U_OPT_OUT_LENGTH   2

/* index names for readability (optional but avoids magic numbers) */
#define IDX_ID      0
#define IDX_IQ      1
#define IDX_OMEGA   2
#define IDX_THETA   3
#define IDX_LD      4
#define IDX_LQ      5
#define IDX_TL      6
#define IDX_WREF    7

/* --- software / C-simulation types (verify against MATLAB first) --- */
typedef double data_t_x_hat_in;
typedef double data_t_u_opt_out;

/* --- HLS fixed-point variant (enable later, after the double version passes) ---
 * #include "ap_fixed.h"
 * typedef ap_fixed<32,16> data_t_x_hat_in;
 * typedef ap_fixed<32,16> data_t_u_opt_out;
 */

#endif /* FOO_DATA_H */
