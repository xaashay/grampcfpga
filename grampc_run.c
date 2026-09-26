/* This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright (C) 2014-2018 by Tobias Englert, Knut Graichen, Felix Mesmer,
 * Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * Developed at the Institute of Measurement, Control, and Microtechnology,
 * Ulm University. All rights reserved.
 *
 * GRAMPC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * GRAMPC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GRAMPC. If not, see <http://www.gnu.org/licenses/>
 *
 */

#ifdef __cplusplus
extern "C" {
#endif


//#include "euler1.h"
#include "grampc_run.h"
/* Heun integrator (grampc_intsys_bridge.c) -- replaces the hardcoded Euler.
 * Euler diverges to NaN for w>=200 here: h/tau_d = 3.0ms/3.37ms = 0.89 is at
 * its stability edge, and the speed-coupling terms push it over.
 * SpeedMPC_Setting.m specifies Integrator='heun' for exactly this reason. */
void intsysHeun_wsys_impl(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t,
                          ctypeRNum *x, ctypeRNum *u, ctypeRNum *p_, typeGRAMPC *grampc);
void intsysHeun_wadjsys_impl(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t,
                             ctypeRNum *x, ctypeRNum *u, ctypeRNum *p_, typeGRAMPC *grampc);

void grampc_run(typeGRAMPC *grampc)
{
	typeInt i, j, k, imult, igrad;
	typeBoolean sysintegrated = 0;
	typeRNum alpha;
	typeRNum cfct_norm;
	typeRNum pen_norm;

	typeRNum *t = grampc->t;
	typeRNum *u = grampc->u;
	typeRNum *gradu = grampc->gradu;
	typeRNum *p = grampc->p;
//	typeRNum *gradp = grampc->gradp;
	typeRNum *T = &grampc->T;
	typeRNum *gradT = &grampc->gradT;

	/* Reset solution structure */
	typeBoolean converged_grad = 0;
	typeBoolean converged_const = 0;
	for (imult = 0; imult < grampc->MaxMultIter; imult++) {
		grampc->iter[imult] = 0;
	}
	//grampc->status = STATUS_NONE;

	/* Validate necessary parameters */
/*	if (grampc->dt <= 0.0) {
		grampc_error(DT_NOT_VALID);
	}
	if (grampc->Thor < grampc->dt) {
		grampc_error(THOR_NOT_VALID);
	}*/

	/* Initial condition */
	if (grampc->ScaleProblem == INT_ON) {
		scale_states(grampc->x, grampc->x0, grampc);
	}
	else {
		MatCopy(grampc->x, grampc->x0, 1, grampc->Nx);
	}

	/* Shift of input, Lagrange multiplier, and penalty parameter trajectories */
	/* Attention shiftTrajecotry and shortenTrajectory require a uniform grid size over the horizon */
	if (grampc->ShiftControl == INT_ON) {
		if (grampc->OptimTime == INT_OFF) {
			shiftTrajectory(u, grampc->Nhor, grampc->Nu, grampc->Nu, grampc->dt, t);
			shiftTrajectory(grampc->uprev, grampc->Nhor, grampc->Nu, grampc->Nu, grampc->dt, t);
			shiftTrajectory(grampc->graduprev, grampc->Nhor, grampc->Nu, grampc->Nu, grampc->dt, t);
			shiftTrajectory(grampc->cfctprev, grampc->Nhor, grampc->Nc, (grampc->Ng + grampc->Nh), grampc->dt, t);
			shiftTrajectory(grampc->mult, grampc->Nhor, grampc->Nc, (grampc->Ng + grampc->Nh), grampc->dt, t);
			shiftTrajectory(grampc->rws_pen, grampc->Nhor, grampc->Nc, (grampc->Ng + grampc->Nh), grampc->dt, t);
		}
		else {
			T[0] = T[0] < grampc->Tmin + grampc->dt ? grampc->Tmin : T[0] - grampc->dt;
			shortenTrajectory(u, grampc->Nhor, grampc->Nu, grampc->Nu, grampc->dt, t);
			shortenTrajectory(grampc->uprev, grampc->Nhor, grampc->Nu, grampc->Nu, grampc->dt, t);
			shortenTrajectory(grampc->graduprev, grampc->Nhor, grampc->Nu, grampc->Nu, grampc->dt, t);
			shortenTrajectory(grampc->cfctprev, grampc->Nhor, grampc->Nc, (grampc->Ng + grampc->Nh), grampc->dt, t);
			shortenTrajectory(grampc->mult, grampc->Nhor, grampc->Nc, (grampc->Ng + grampc->Nh), grampc->dt, t);
			shortenTrajectory(grampc->rws_pen, grampc->Nhor, grampc->Nc, (grampc->Ng + grampc->Nh), grampc->dt, t);
			//discretize_time(t, T[0], &grampc);

			typeInt i;
            typeInt N = (grampc)->Nhor;
            typeRNum m, c;
            typeRNum T_ = T[0];
            typeRNum *tvec = t;

			if ((grampc)->TimeDiscretization == INT_UNIFORM || T_ <= (N - 1) * (grampc)->dt) {
                    m = T_ / (N - 1);
                    for (i = 0; i < (grampc)->Nhor; i++) {
                        tvec[i] = m * i;
                }
            }
	/* Nonuniform discretization of interval [0, T_] */
            else if ((grampc)->TimeDiscretization == INT_NONUNIFORM) {
                m = (T_ / (N - 1) - (grampc)->dt) / (N - 2);
                c = (grampc)->dt - m;
                for (i = 0; i < (grampc)->Nhor; i++) {
                    tvec[i] = m * i * i + c * i;
                }
            }

			}
		}


	/* LOOP OVER NO. OF MULTIPLIER STEPS *******************************************/
	for (imult = 0; imult < grampc->MaxMultIter; imult++) {

		/* LOOP OVER NO. OF GRADIENT STEPS *********************************************/
		for (igrad = 0; igrad < grampc->MaxGradIter; igrad++) {

			/* Forward integration of system and evaluation of constraints */
			if (!sysintegrated) {
				evaluate_sys(t, u, p, grampc);
				evaluate_constraints(t, u, p, 1, 0, grampc);
				sysintegrated = 1;
			}

			/* Backward integration of adjoint system */
			evaluate_adjsys(t, u, p, grampc);

			if (grampc->OptimControl == INT_ON) {
				/* Gradient w.r.t. u */
				evaluate_gradu(grampc);
			}
			if (grampc->OptimParam == INT_ON) {
				/* Gradient w.r.t. p */
			//	evaluate_gradp(grampc);
			}
	/*		if (grampc->OptimTime == INT_ON) {
				// Gradient w.r.t. T
				evaluate_gradT(grampc);
			}*/

			/* Determine step size alpha by linesearch */
			alpha = 0;
			if (grampc->LineSearchType == INT_ADAPTIVELS) {
				/* Adaptive line search */
	//			linesearch_adaptive(&alpha, igrad, grampc);
			}
			else {
				/* Explicit line search */
				linesearch_explicit(&alpha, grampc);
			}

			if (grampc->OptimControl == INT_ON) {
				// Save previous values
				MatCopy(grampc->uprev, u, grampc->Nhor, grampc->Nu);
				MatCopy(grampc->graduprev, gradu, grampc->Nhor, grampc->Nu);
				// Update control
				for (i = 0; i < grampc->Nhor; i++) {
					for (j = 0; j < grampc->Nu; j++) {
						u[i*grampc->Nu + j] = u[i*grampc->Nu + j] - alpha * gradu[i*grampc->Nu + j];
					}
				}
				inputproj(u, grampc);
			}
/*			if (grampc->OptimParam == INT_ON) {
				// Save previous values
				MatCopy(grampc->pprev, p, 1, grampc->Np);
				MatCopy(grampc->gradpprev, gradp, 1, grampc->Np);
				// Update parameters
				for (j = 0; j < grampc->Np; j++) {
					p[j] = p[j] - grampc->OptimParamLineSearchFactor * alpha * gradp[j];
				}
				paramproj(p, grampc);
			}*/
			if (grampc->OptimTime == INT_ON) {
				/* Save previous values */
				grampc->Tprev = T[0];
				grampc->gradTprev = gradT[0];
				/* Update time */
				T[0] = T[0] - grampc->OptimTimeLineSearchFactor * alpha * gradT[0];
				timeproj(T, grampc);
			//	discretize_time(t, T[0],&(*grampc));
			typeInt i;
            typeInt N = (grampc)->Nhor;
            typeRNum m, c;
            typeRNum T_ = T[0];
            typeRNum *tvec = t;

			if ((grampc)->TimeDiscretization == INT_UNIFORM || T_ <= (N - 1) * (grampc)->dt) {
                    m = T_ / (N - 1);
                    for (i = 0; i < (grampc)->Nhor; i++) {
                        tvec[i] = m * i;
                }
            }
	/* Nonuniform discretization of interval [0, T_] */
            else if ((grampc)->TimeDiscretization == INT_NONUNIFORM) {
                m = (T_ / (N - 1) - (grampc)->dt) / (N - 2);
                c = (grampc)->dt - m;
                for (i = 0; i < (grampc)->Nhor; i++) {
                    tvec[i] = m * i * i + c * i;
                }
            }

			}

			/* Trajectories updated, integration necessary */
			sysintegrated = 0;

			/* Convergence test for gradient */
			if (grampc->ConvergenceCheck == INT_ON) {
				converged_grad = convergence_test_gradient(grampc->ConvergenceGradientRelTol, grampc);
				if (converged_grad) {
					//grampc->status |= STATUS_GRADIENT_CONVERGED;
					igrad++;
					break;
				}
			}
		}
		/* END GRADIENT LOOP ***********************************************************/
		grampc->iter[imult] = igrad;

		/* Forward integration of system and evaluation of constraints */
		evaluate_sys(t, u, p, grampc);
		evaluate_constraints(t, u, p, (imult + 1 < grampc->MaxMultIter), 1, grampc);
		sysintegrated = 1;

		/* Convergence test for constraints */
		if (grampc->ConvergenceCheck == INT_ON && converged_grad) {
			converged_const = convergence_test_constraints(grampc->cfctAbsTol, grampc);
			if (converged_const) {
				//grampc->status |= STATUS_CONSTRAINTS_CONVERGED;
				break;
			}
		}
	}
	/* END MULTIPLIER LOOP ***********************************************************/

	/* Calculation of xnext */
	interplin(grampc->xnext, t, grampc->x, grampc->dt, grampc->Nx, grampc->Nhor, 1);

	/* Unscaling */
	if (grampc->ScaleProblem == INT_ON) {
		unscale_states(grampc->xnext, grampc->xnext, grampc);
		unscale_controls(grampc->unext, u, grampc);
//		unscale_parameters(grampc->pnext, p, grampc);
		unscale_time(&grampc->Tnext, T[0], grampc);
	}
	else {
		MatCopy(grampc->unext, u, 1, grampc->Nu);
		//MatCopy(grampc->pnext, p, 1, grampc->Np);
		grampc->Tnext = T[0];
	}

	/* Calculation of cost J */
	evaluate_cost(grampc->J, t, u, p, grampc);
	/* Calculation of constraint norm and penalty norm*/
	cfct_norm = 0;
	for (i = 0; i < grampc->Nhor; i++) {
		k = i * grampc->Nc;
		for (j = 0; j < grampc->Ng; j++) {
			cfct_norm = cfct_norm + grampc->rws_cfct[k + j] * grampc->rws_cfct[k + j];
		}
		for (; j < grampc->Ng + grampc->Nh; j++) {
			cfct_norm = cfct_norm + MAX(grampc->rws_cfct[k + j], 0) * MAX(grampc->rws_cfct[k + j], 0);
		}
		for (; j < grampc->Ng + grampc->Nh + grampc->NgT; j++) {
			cfct_norm = cfct_norm + grampc->rws_cfct[k + j] * grampc->rws_cfct[k + j];
		}
		for (; j < grampc->Ng + grampc->Nh + grampc->NgT + grampc->NhT; j++) {
			cfct_norm = cfct_norm + MAX(grampc->rws_cfct[k + j], 0) * MAX(grampc->rws_cfct[k + j], 0);
		}
	}
	cfct_norm = SQRT(cfct_norm);
	MatNorm(&pen_norm, grampc->rws_pen, grampc->Nhor, grampc->Nc);
	if (cfct_norm > grampc->cfct && pen_norm >= grampc->pen && !converged_const) {
		//grampc->status |= STATUS_INFEASIBLE;
	}
	grampc->cfct = cfct_norm;
	grampc->pen = pen_norm;
}


void evaluate_constraints(ctypeRNum *t, ctypeRNum *u, ctypeRNum *p, const typeBoolean evaljac, const typeBoolean updatemultiplier, typeGRAMPC *grampc)
{
	typeInt i;
    typeInt N;

	ctypeRNum *x_ = NULL;
	ctypeRNum *u_ = NULL;
	ctypeRNum *p_ = p;
	typeRNum *mult = NULL;
	typeRNum *pen = NULL;
	typeRNum *cfct = NULL;
	typeRNum *cfctprev = NULL;
	typeRNum *dcdx = NULL;
	typeRNum *dcdu = NULL;
	//typeRNum *dcdp = NULL;
	typeRNum *dcdt = &grampc->dcdt;
	typeRNum *thresholds = NULL;
	typeRNum *cScale = NULL;
	typeBoolean converged_grad = 0;

	typeRNum *c = grampc->rwsGeneral; /* size:  Nc+ 2*(Nx+Nu+Np) */
	typeRNum *dgdxvec = c + grampc->Nc;
	typeRNum *dhdxvec = dgdxvec + grampc->Nx;
	typeRNum *dgdpvec = dhdxvec + grampc->Nx;
	typeRNum *dhdpvec = dgdpvec + grampc->Np;
	typeRNum *dgduvec = dhdpvec + grampc->Np;
	typeRNum *dhduvec = dgduvec + grampc->Nu;
	typeRNum dgTdT, dhTdT;

	/* return if no constraints are defined */
	if (grampc->Nc == 0) {
		return;
	}

	/* Init constraint jacobians */
	MatSetScalar(grampc->dcdx, 0, grampc->Nhor + 1, grampc->Nx);
	MatSetScalar(grampc->dcdu, 0, grampc->Nhor, grampc->Nu);
//	MatSetScalar(grampc->dcdp, 0, grampc->Nhor + 1, grampc->Np);
	MatSetScalar(dgdxvec, 0, 1, 2 * (grampc->Nx + grampc->Np + grampc->Nu));

	/* Multiplier update depends on the convergence of the subproblem */
	if (updatemultiplier) {
		converged_grad = convergence_test_gradient(grampc->AugLagUpdateGradientRelTol, grampc);
		if (converged_grad) {
			//grampc->status |= STATUS_MULTIPLIER_UPDATE;
		}
	}

	/* Unscaling */
	if (grampc->ScaleProblem == INT_ON) {
		ASSIGN_P(p_, p, grampc);
	}

	/* loop over the prediction horizon */
	if ((grampc->Ng + grampc->Nh > 0) && ((grampc->EqualityConstraints == INT_ON) || (grampc->InequalityConstraints == INT_ON))) {
        /* if there are terminal constraints, the integral constraints are not evaluated for the last point */
        if (grampc->NgT + grampc->NhT > 0) {
            N = grampc->Nhor - 1;
        }
        /* if there are no terminal constraints, the integral constraints are evaluated for all points */
        else {
            N = grampc->Nhor;
        }
        for (i = 0; i < N; i++)
		{
			mult = grampc->mult + i * grampc->Nc;
			pen = grampc->rws_pen + i * grampc->Nc;
			cfct = grampc->rws_cfct + i * grampc->Nc;
			cfctprev = grampc->cfctprev + i * grampc->Nc;
			dcdx = grampc->dcdx + i * grampc->Nx;
			dcdu = grampc->dcdu + i * grampc->Nu;
//			dcdp = grampc->dcdp + i * grampc->Np;
			thresholds = grampc->cfctAbsTol;
			cScale = grampc->cScale;

			/* Unscaling */
			if (grampc->ScaleProblem == INT_ON) {
				ASSIGN_X(x_, grampc->x + i * grampc->Nx, grampc);
				ASSIGN_U(u_, u + i * grampc->Nu, grampc);
			}
			else {
				x_ = grampc->x + i * grampc->Nx;
				u_ = u + i * grampc->Nu;
			}

			/* EQUALITY CONSTRAINTS **************************************************/
			if (grampc->Ng > 0 && grampc->EqualityConstraints == INT_ON) {

				/* Evaluate constraints */
				gfct(cfct, t[i], x_, u_, p_, grampc->userparam);
				if (grampc->ScaleProblem == INT_ON) {
					scale_constraints(cfct, cScale, grampc->Ng);
				}

				/* Update multipliers */
				if (updatemultiplier) {
					update_multiplier_eqc(mult, pen, cfct, cfctprev, thresholds, grampc->Ng, converged_grad, grampc);
				}

				/* Evaluate jacobians */
				if (evaljac) {
					compute_jacobian_multiplier(c, mult, pen, cfct, grampc->Ng);
					if (grampc->ScaleProblem == INT_ON) {
						scale_constraints(c, cScale, grampc->Ng);
					}

					dgdx_vec(dgdxvec, t[i], x_, u_, p_, c, grampc->userparam);
					MatAdd(dcdx, dcdx, dgdxvec, 1, grampc->Nx);

					if (grampc->OptimControl == INT_ON) {
						dgdu_vec(dgduvec, t[i], x_, u_, p_, c, grampc->userparam);
						MatAdd(dcdu, dcdu, dgduvec, 1, grampc->Nu);
					}

			/*		if (grampc->OptimParam == INT_ON) {
						dgdp_vec(dgdpvec, t[i], x_, u_, p_, c, grampc->userparam);
						MatAdd(dcdp, dcdp, dgdpvec, 1, grampc->Np);
					}*/
				}
			}

			mult = mult + grampc->Ng;
			pen = pen + grampc->Ng;
			cfct = cfct + grampc->Ng;
			cfctprev = cfctprev + grampc->Ng;
			thresholds = thresholds + grampc->Ng;
			cScale = cScale + grampc->Ng;

			/* INEQUALITY CONSTRAINTS ************************************************/
			if (grampc->Nh > 0 && grampc->InequalityConstraints == INT_ON) {

				/* Evaluate constraints */
				hfct(cfct, t[i], x_, u_, p_, grampc->userparam);
				if (grampc->ScaleProblem == INT_ON) {
					scale_constraints(cfct, cScale, grampc->Nh);
				}
				update_cfct_for_ieqc(mult, pen, cfct, grampc->Nh);

				/* Update multipliers */
				if (updatemultiplier) {
					update_multiplier_ieqc(mult, pen, cfct, cfctprev, thresholds, grampc->Nh, converged_grad, grampc);
				}

				/* Evaluate jacobians */
				if (evaljac) {
					compute_jacobian_multiplier(c, mult, pen, cfct, grampc->Nh);
					if (grampc->ScaleProblem == INT_ON) {
						scale_constraints(c, cScale, grampc->Nh);
					}

					dhdx_vec(dhdxvec, t[i], x_, u_, p_, c, grampc->userparam);
					MatAdd(dcdx, dcdx, dhdxvec, 1, grampc->Nx);

					if (grampc->OptimControl == INT_ON) {
						dhdu_vec(dhduvec, t[i], x_, u_, p_, c, grampc->userparam);
						MatAdd(dcdu, dcdu, dhduvec, 1, grampc->Nu);
					}

/*					if (grampc->OptimParam == INT_ON) {
						dhdp_vec(dhdpvec, t[i], x_, u_, p_, c, grampc->userparam);
						MatAdd(dcdp, dcdp, dhdpvec, 1, grampc->Np);
					}*/
				}
			}
		}
	}

	/* Terminal Constraints */
	if ((grampc->NgT + grampc->NhT > 0) && ((grampc->TerminalEqualityConstraints == INT_ON) || (grampc->TerminalInequalityConstraints == INT_ON)) ){

		i = grampc->Nhor - 1;
		mult = grampc->mult + i * grampc->Nc + grampc->Ng + grampc->Nh;
		pen = grampc->rws_pen + i * grampc->Nc + grampc->Ng + grampc->Nh;
		cfct = grampc->rws_cfct + i * grampc->Nc + grampc->Ng + grampc->Nh;
		cfctprev = grampc->cfctprev + i * grampc->Nc + grampc->Ng + grampc->Nh;
		/* dcdx for Terminal constraint must be different from dcdx for integral constraints */
		dcdx = grampc->dcdx + grampc->Nhor * grampc->Nx;
		/* dcdp for Terminal constraint must be different from dcdp for integral constraints */
//		dcdp = grampc->dcdp + grampc->Nhor * grampc->Np;
		/* dcdt for Terminal constraint must be different from dcdt for integral constraints */
		*dcdt = 0;
		thresholds = grampc->cfctAbsTol + grampc->Ng + grampc->Nh;
		cScale = grampc->cScale + grampc->Ng + grampc->Nh;

		/* Init constraint jacobians */;
		MatSetScalar(dgdxvec, 0, 1, 2 * (grampc->Nx + grampc->Np));
		dgTdT = 0;
		dhTdT = 0;

		/* Unscaling */
		if (grampc->ScaleProblem == INT_ON) {
			ASSIGN_X(x_, grampc->x + i * grampc->Nx, grampc);
		}
		else {
			x_ = grampc->x + i * grampc->Nx;
		}

		/* TERMINAL EQUALITY CONSTRAINTS *******************************************/
		if (grampc->NgT > 0 && grampc->TerminalEqualityConstraints == INT_ON) {

			/* Evaluate constraints */
			gTfct(cfct, t[i], x_, p_, grampc->userparam);
			if (grampc->ScaleProblem == INT_ON) {
				scale_constraints(cfct, cScale, grampc->NgT);
			}

			/* Update multipliers */
			if (updatemultiplier) {
				update_multiplier_eqc(mult, pen, cfct, cfctprev, thresholds, grampc->NgT, converged_grad, grampc);
			}

			/* Evaluate jacobians */
			if (evaljac) {
				compute_jacobian_multiplier(c, mult, pen, cfct, grampc->NgT);
				if (grampc->ScaleProblem == INT_ON) {
					scale_constraints(c, cScale, grampc->NgT);
				}

				dgTdx_vec(dgdxvec, t[i], x_, p_, c, grampc->userparam);
				MatAdd(dcdx, dcdx, dgdxvec, 1, grampc->Nx);

/*				if (grampc->OptimParam == INT_ON) {
					dgTdp_vec(dgdpvec, t[i], x_, p_, c, grampc->userparam);
					MatAdd(dcdp, dcdp, dgdpvec, 1, grampc->Np);
				}*/

				if (grampc->OptimTime == INT_ON) {
					dgTdT_vec(&dgTdT, t[i], x_, p_, c, grampc->userparam);
					*dcdt = *dcdt + dgTdT;
				}
			}
		}

		mult = mult + grampc->NgT;
		pen = pen + grampc->NgT;
		cfct = cfct + grampc->NgT;
		cfctprev = cfctprev + grampc->NgT;
		thresholds = thresholds + grampc->NgT;
		cScale = cScale + grampc->NgT;

		/* TERMINAL INEQUALITY CONSTRAINTS *****************************************/
		if (grampc->NhT > 0 && grampc->TerminalInequalityConstraints == INT_ON) {

			/* Evaluate constraints */
			hTfct(cfct, t[i], x_, p_, grampc->userparam);
			if (grampc->ScaleProblem == INT_ON) {
				scale_constraints(cfct, cScale, grampc->NhT);
			}
			update_cfct_for_ieqc(mult, pen, cfct, grampc->NhT);

			/* Update multipliers */
			if (updatemultiplier) {
				update_multiplier_ieqc(mult, pen, cfct, cfctprev, thresholds, grampc->NhT, converged_grad, grampc);
			}

			/* Evaluate jacobians */
			if (evaljac) {
				compute_jacobian_multiplier(c, mult, pen, cfct, grampc->NhT);
				if (grampc->ScaleProblem == INT_ON) {
					scale_constraints(c, cScale, grampc->NhT);
				}

				dhTdx_vec(dhdxvec, t[i], x_, p_, c, grampc->userparam);
				MatAdd(dcdx, dcdx, dhdxvec, 1, grampc->Nx);

/*				if (grampc->OptimParam == INT_ON) {
					dhTdp_vec(dhdpvec, t[i], x_, p_, c, grampc->userparam);
					MatAdd(dcdp, dcdp, dhdpvec, 1, grampc->Np);
				}*/

				if (grampc->OptimTime == INT_ON) {
					dhTdT_vec(&dhTdT, t[i], x_, p_, c, grampc->userparam);
					*dcdt = *dcdt + dhTdT;
				}
			}
		}
	}
}

void update_multiplier_eqc(typeRNum *mult, typeRNum *pen, ctypeRNum *cfct, typeRNum *cfctprev,
	ctypeRNum *thresholds, ctypeInt Ncon, typeBoolean converged_grad, typeGRAMPC *grampc)
{
	typeInt i;
	for (i = 0; i < Ncon; i++) {
		/* Increase multipliers for violated constraints if minimization converged */
		if (ABS(cfct[i]) > thresholds[i] && converged_grad) {
			if (grampc->ConstraintsHandling == INT_AUGLAG) {
				mult[i] = mult[i] + (1 - grampc->MultiplierDampingFactor) * pen[i] * cfct[i];
				if (mult[i] > grampc->MultiplierMax) {
					mult[i] = grampc->MultiplierMax;
					//grampc->status |= STATUS_MULTIPLIER_MAX;
				}
				else if (mult[i] < -grampc->MultiplierMax) {
					mult[i] = -grampc->MultiplierMax;
					//grampc->status |= STATUS_MULTIPLIER_MAX;
				}
			}
			if (ABS(cfct[i]) > grampc->PenaltyIncreaseThreshold * ABS(cfctprev[i])) {
				pen[i] = pen[i] * grampc->PenaltyIncreaseFactor;
				if (pen[i] > grampc->PenaltyMax) {
					pen[i] = grampc->PenaltyMax;
					//grampc->status |= STATUS_PENALTY_MAX;
				}
			}
			/* Save cfct as cfctprev for next multiplier update */
			cfctprev[i] = cfct[i];
		}
		/* Decrease multipliers for satisfied constraints */
		if (ABS(cfct[i]) < thresholds[i] / 10.0) {
			pen[i] = MAX(pen[i] * grampc->PenaltyDecreaseFactor, grampc->PenaltyMin);
		}
	}
}

void update_multiplier_ieqc(typeRNum *mult, typeRNum *pen, ctypeRNum *cfct, typeRNum *cfctprev,
	ctypeRNum *thresholds, ctypeInt Ncon, typeBoolean converged_grad, typeGRAMPC *grampc)
{
	typeInt i;
	for (i = 0; i < Ncon; i++) {
		/* Increase multipliers for violated constraints if minimization converged */
		if (cfct[i] > thresholds[i] && converged_grad) {
			if (grampc->ConstraintsHandling == INT_AUGLAG) {
				mult[i] = mult[i] + (1 - grampc->MultiplierDampingFactor) * pen[i] * cfct[i];
				if (mult[i] > grampc->MultiplierMax) {
					mult[i] = grampc->MultiplierMax;
					//grampc->status |= STATUS_MULTIPLIER_MAX;
				}
			}
			if (cfct[i] > grampc->PenaltyIncreaseThreshold * cfctprev[i]) {
				pen[i] = pen[i] * grampc->PenaltyIncreaseFactor;
				if (pen[i] > grampc->PenaltyMax) {
					pen[i] = grampc->PenaltyMax;
					//grampc->status |= STATUS_PENALTY_MAX;
				}
			}
			/* Save cfct as cfctprev for next multiplier update */
			cfctprev[i] = cfct[i];
		}
		/* Decrease multipliers for satisfied constraints */
		if (cfct[i] < thresholds[i] / 10.0) {
			if (grampc->ConstraintsHandling == INT_AUGLAG && cfct[i] < 0) {
				mult[i] = mult[i] + (1 - grampc->MultiplierDampingFactor) * pen[i] * cfct[i];
			}
			pen[i] = MAX(pen[i] * grampc->PenaltyDecreaseFactor, grampc->PenaltyMin);
		}
	}
}

void update_cfct_for_ieqc(ctypeRNum *mult, ctypeRNum *pen, typeRNum *cfct, ctypeInt Ncon)
{
	typeInt i;
	for (i = 0; i < Ncon; i++) {
		cfct[i] = MAX(cfct[i], -mult[i] / pen[i]);
	}
}

void compute_jacobian_multiplier(typeRNum *c, ctypeRNum *mult, ctypeRNum *pen, ctypeRNum *cfct, ctypeInt Ncon)
{
	typeInt i;
	for (i = 0; i < Ncon; i++) {
		c[i] = mult[i] + pen[i] * cfct[i];
	}
}


void evaluate_sys(ctypeRNum *t, ctypeRNum *u, ctypeRNum *p, typeGRAMPC *grampc)
{
//	typeIntffctPtr pIntSys;
	ctypeRNum *p_ = p;

	/* integrator for system integration */
/*	if (grampc->Integrator == INT_EULER) {
		//pIntSys = &intsysEuler_wsys;
		intsysEuler_wsys(grampc->x, FWINT, grampc->Nhor, t, grampc->x, u, p_, grampc);
	}
	else if (grampc->Integrator == INT_MODEULER) {
		pIntSys = &intsysModEuler;
	}
	else if (grampc->Integrator == INT_HEUN) {
		pIntSys = &intsysHeun;
	}
	else if (grampc->Integrator == INT_RODAS) {
		pIntSys = &intsysRodas;
	}
	else {
		pIntSys = &intsysRuKu45;
	}
*/
	/* Unscaling */
	if (grampc->ScaleProblem == INT_ON) {
		ASSIGN_P(p_, p, grampc);
	}

	//(*pIntSys)(grampc->x, FWINT, grampc->Nhor, t, grampc->x, u, p_, grampc, &Wsys);
	intsysHeun_wsys_impl(grampc->x, FWINT, grampc->Nhor, t, grampc->x, u, p_, grampc);
}

void evaluate_adjsys(ctypeRNum *t, ctypeRNum *u, ctypeRNum *p, typeGRAMPC *grampc)
{
	typeInt i = grampc->Nhor - 1;
	ctypeRNum *x_ = grampc->x + i * grampc->Nx;
	ctypeRNum *p_ = p;
	typeRNum *adj_ = grampc->adj + i * grampc->Nx;

	//typeIntffctPtr pIntSys;

	/* integrator for adjoint system integration */
/*	if (grampc->Integrator == INT_EULER) {
		//pIntSys = &intsysEuler_wadjsys;
		intsysEuler_wadjsys(adj_, BWINT, grampc->Nhor, t + i, grampc->x + i * grampc->Nx, u + i * grampc->Nu, p_, grampc);
	}
	else if (grampc->Integrator == INT_MODEULER) {
		pIntSys = &intsysModEuler;
	}
	else if (grampc->Integrator == INT_HEUN) {
		pIntSys = &intsysHeun;
	}
	else if (grampc->Integrator == INT_RODAS) {
		pIntSys = &intsysRodas;
	}
	else {
		pIntSys = &intsysRuKu45;
	}
*/
	/* Unscaling */
	if (grampc->ScaleProblem == INT_ON) {
		ASSIGN_X(x_, grampc->x + i * grampc->Nx, grampc);
		ASSIGN_P(p_, p, grampc);
	}

	/* Terminal condition for adjoint states */
	MatSetScalar(adj_, 0, 1, grampc->Nx);

	/* Jacobian of Terminal cost */
	if (grampc->TerminalCost == INT_ON) {
		dVdx(adj_, t[i], x_, p_, grampc->xdes, grampc->userparam);
/*		if (grampc->ScaleProblem == INT_ON) {
			scale_cost(adj_, grampc->JScale, grampc->Nx);
		}*/
	}

	/* Jacobian of Terminal equality and inequality constraints */
	if ((grampc->NgT + grampc->NhT > 0) && (grampc->TerminalEqualityConstraints == INT_ON
		|| grampc->TerminalInequalityConstraints == INT_ON)) {
		MatAdd(adj_, adj_, grampc->dcdx + grampc->Nhor * grampc->Nx, 1, grampc->Nx);
	}

	/* Scaling */
	if (grampc->ScaleProblem == INT_ON) {
		scale_adjoints(adj_, adj_, grampc);
	}

	/* integration of adjoint system in reverse time */
	//(*pIntSys)(adj_, BWINT, grampc->Nhor, t + i, grampc->x + i * grampc->Nx, u + i * grampc->Nu, p_, grampc, &Wadjsys);
	intsysHeun_wadjsys_impl(adj_, BWINT, grampc->Nhor, t + i, grampc->x + i * grampc->Nx, u + i * grampc->Nu, p_, grampc);
}

void Wsys(typeRNum *s, ctypeRNum *x, ctypeRNum *t, ctypeRNum *dummy,
	ctypeRNum *u, ctypeRNum *p_, ctypeRNum *dcdx, typeGRAMPC *grampc)
{
	ctypeRNum *x_ = x;
	ctypeRNum *u_ = u;
	typeInt i;

	/* remove warning unreferenced formal parameter */
	(void)(dummy);
	(void)(dcdx);

	/* Unscaling */
	if (grampc->ScaleProblem == INT_ON) {
		ASSIGN_X(x_, x, grampc);
		ASSIGN_U(u_, u, grampc);
	}

	//ffct(s, t[0] + grampc->t0, x_, u_, p_, grampc->userparam);
	ffct(s, t[0] + grampc->t0, x_, u_, grampc->userparam);

	/* Scaling */
	if (grampc->ScaleProblem == INT_ON) {
		for (i = 0; i < grampc->Nx; i++) {
			s[i] = s[i] / grampc->xScale[i];
		}
	}
}

void Wadjsys(typeRNum *s, ctypeRNum *adj, ctypeRNum *t, ctypeRNum *x,
	ctypeRNum *u, ctypeRNum *p_, ctypeRNum *dcdx, typeGRAMPC *grampc)
{
	ctypeRNum *x_ = x;
	ctypeRNum *adj_ = adj;
	ctypeRNum *u_ = u;
	typeRNum *dLdx = grampc->rwsGeneral; /* size:  Nx */

	typeInt i;

	/* Unscaling */
	if (grampc->ScaleProblem == INT_ON) {
		ASSIGN_X(x_, x, grampc);
		ASSIGN_ADJ(adj_, adj, grampc);
		ASSIGN_U(u_, u, grampc);
	}

	MatSetScalar(dLdx, 0, 1, grampc->Nx);

	if (grampc->IntegralCost == INT_ON) {
		dldx(dLdx, t[0], x_, u_, p_, grampc->xdes, grampc->udes, grampc->userparam);
/*		if (grampc->ScaleProblem == INT_ON) {
			scale_cost(dLdx, grampc->JScale, grampc->Nx);
		}*/
	}

	dfdx_vec(s, t[0] + grampc->t0, x_, adj_, u_, p_, grampc->userparam);

	/* Scaling */
	if (grampc->ScaleProblem == INT_ON) {
		for (i = 0; i < grampc->Nx; i++) {
			s[i] = (-dLdx[i] - s[i] - dcdx[i]) * grampc->xScale[i];
		}
	}
	else {
		for (i = 0; i < grampc->Nx; i++) {
			s[i] = -dLdx[i] - s[i] - dcdx[i];
		}
	}
}

void evaluate_gradu(typeGRAMPC *grampc)
{
	typeInt i, j;

	ctypeRNum *t = grampc->t;
	ctypeRNum *x_ = NULL;
	ctypeRNum *u_ = NULL;
	ctypeRNum *p_ = grampc->p;
	ctypeRNum *adj_ = NULL;
	ctypeRNum *dcdu = NULL;

	typeRNum *dLdu = grampc->rwsGeneral;
	typeRNum *s = dLdu + grampc->Nu;
	/*                      dLdu  s  */
	/* sizeof(rwsGeneral) = Nu + Nu */

	MatSetScalar(dLdu, 0, 1, grampc->Nu);

	/* Unscaling */
/*	if (grampc->ScaleProblem == INT_ON) {
		ASSIGN_P(p_, grampc->p, grampc);
	}*/

	for (i = 0; i < grampc->Nhor; i++) {

		dcdu = grampc->dcdu + i * grampc->Nu;

		/* Unscaling */
		if (grampc->ScaleProblem == INT_ON) {
			ASSIGN_X(x_, grampc->x + i * grampc->Nx, grampc);
			ASSIGN_ADJ(adj_, grampc->adj + i * grampc->Nx, grampc);
			ASSIGN_U(u_, grampc->u + i * grampc->Nu, grampc);
		}
		else {
			x_ = grampc->x + i * grampc->Nx;
			adj_ = grampc->adj + i * grampc->Nx;
			u_ = grampc->u + i * grampc->Nu;
		}

		if (grampc->IntegralCost == INT_ON) {
			dldu(dLdu, t[i], x_, u_, p_, grampc->xdes, grampc->udes, grampc->userparam);
	/*		if (grampc->ScaleProblem == INT_ON) {
				scale_cost(dLdu, grampc->JScale, grampc->Nu);
			}*/
		}

		dfdu_vec(s, t[i] + grampc->t0, x_, adj_, u_, p_, grampc->userparam);

		/* Scaling */
		if (grampc->ScaleProblem == INT_ON) {
			for (j = 0; j < grampc->Nu; j++) {
				grampc->gradu[i*grampc->Nu + j] = (dLdu[j] + s[j] + dcdu[j]) * grampc->uScale[j];
			}
		}
		else {
			for (j = 0; j < grampc->Nu; j++) {
				grampc->gradu[i*grampc->Nu + j] = dLdu[j] + s[j] + dcdu[j];
			}
		}
	}
}

void inputproj(typeRNum *u, typeGRAMPC *grampc)
{
	typeInt i, j;
	typeRNum *umin_ = grampc->umin;
	typeRNum *umax_ = grampc->umax;

	if (grampc->ScaleProblem == INT_ON) {
		umin_ = grampc->rwsScale + 2 * grampc->Nx;
		umax_ = umin_ + grampc->Nu;
		scale_controls(umin_, grampc->umin, grampc);
		scale_controls(umax_, grampc->umax, grampc);
	}

	for (i = 0; i < grampc->Nhor; i++) {
		for (j = 0; j < grampc->Nu; j++) {
			/* lower bound */
			if (u[j + i * grampc->Nu] < umin_[j]) {
				u[j + i * grampc->Nu] = umin_[j];
			}
			/* upper bound */
			else if (u[j + i * grampc->Nu] > umax_[j]) {
				u[j + i * grampc->Nu] = umax_[j];
			}
		}
	}
}

/*
void evaluate_gradp(typeGRAMPC *grampc)
{
	typeInt i, j;
	typeRNum h;

	ctypeRNum *t = grampc->t;
	typeRNum *x = NULL;
	typeRNum *adj = NULL;
	typeRNum *u = NULL;
    typeRNum *p_ = grampc->p;
	ctypeRNum *dcdp = NULL;

	//typeRNum *gradp = grampc->gradp;
	typeRNum *gradp = NULL;
	typeRNum *s = grampc->rwsGeneral;

	// Initialize to zero
//	MatSetScalarp(gradp, 0.0, 1, grampc->Np);

	// Unscaling
 	if (grampc->ScaleProblem == INT_ON) {
		ASSIGN_P(p_, grampc->p, grampc);
	}

	// Compute integral over gradp using trapezoidal rule
	for (i = 0; i < grampc->Nhor; i++)
	{
		// Integration step size
		if (i == 0) {
			h = (t[i + 1] - t[i]) / 2;
		}
		else if (i <= grampc->Nhor - 2) {
			h = (t[i + 1] - t[i - 1]) / 2;
		}
		else {
			h = (t[i] - t[i - 1]) / 2;
		}

		x = grampc->x + i * grampc->Nx;
		adj = grampc->adj + i * grampc->Nx;
		u = grampc->u + i * grampc->Nu;
		dcdp = grampc->dcdp + i * grampc->Np;

		WintParam(s, t[i], x, adj, u, p_, dcdp, grampc);
		for (j = 0; j < grampc->Np; j++) {
			gradp[j] = gradp[j] + h * s[j];
		}
	}

	x = grampc->x + (grampc->Nhor - 1) * grampc->Nx;
	dcdp = grampc->dcdp + grampc->Nhor * grampc->Np;

	WtermParam(s, t[i], x, p_, dcdp, grampc);
	MatAdd(gradp, gradp, s, 1, grampc->Np);

	// Scaling
 	if (grampc->ScaleProblem == INT_ON) {
		for (j = 0; j < grampc->Np; j++) {
			gradp[j] = gradp[j] * grampc->pScale[j];
		}
	}
}
*/
void WintParam(typeRNum *s, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj,
	ctypeRNum *u, ctypeRNum *p_, ctypeRNum *dcdp, typeGRAMPC *grampc)
{
	typeInt i;
	ctypeRNum *x_ = x;
	ctypeRNum *adj_ = adj;
	ctypeRNum *u_ = u;

	typeRNum *dldp_val = grampc->rwsGeneral + grampc->Np;
	typeRNum *dfdp_val = dldp_val + grampc->Np;
	MatSetScalar(dldp_val, 0, 1, grampc->Np);

	/* Unscaling */
	if (grampc->ScaleProblem == INT_ON) {
		ASSIGN_X(x_, x, grampc);
		ASSIGN_ADJ(adj_, adj, grampc);
		ASSIGN_U(u_, u, grampc);
	}

	/* dl(x,u,p)/dp + df(x,u,p)/dp' * adj + dc(x,u,p)/dp' * mult */
	if (grampc->IntegralCost == INT_ON) {
		dldp(dldp_val, t, x_, u_, p_, grampc->xdes, grampc->udes, grampc->userparam);
/*		if (grampc->ScaleProblem == INT_ON) {
			scale_cost(dldp_val, grampc->JScale, grampc->Np);
		}*/
	}

	dfdp_vec(dfdp_val, t + grampc->t0, x_, adj_, u_, p_, grampc->userparam);

	for (i = 0; i < grampc->Np; i++) {
		s[i] = dldp_val[i] + dfdp_val[i] + dcdp[i];
	}
}

void WtermParam(typeRNum *s, ctypeRNum t, ctypeRNum *x, ctypeRNum *p_, ctypeRNum *dcdp, typeGRAMPC *grampc)
{
	typeInt i;
	ctypeRNum *x_ = x;

	typeRNum *dldp = grampc->rwsGeneral + grampc->Np;
	MatSetScalar(dldp, 0, 1, grampc->Np);

	/* Unscaling */
	if (grampc->ScaleProblem == INT_ON) {
		ASSIGN_X(x_, x, grampc);
	}

	/* dV(x(T),p)/dp + dc(x(T),p)/dp' * mult(T) */
	if (grampc->TerminalCost == INT_ON) {
		dVdp(dldp, t, x_, p_, grampc->xdes, grampc->userparam);
/*		if (grampc->ScaleProblem == INT_ON) {
			scale_cost(dldp, grampc->JScale, grampc->Np);
		}   */
	}

	for (i = 0; i < grampc->Np; i++) {
		s[i] = dldp[i] + dcdp[i];
	}
}
/*
void paramproj(typeRNum *p, typeGRAMPC *grampc)
{
	typeInt i;
	typeRNum *pmin_ = grampc->pmin;
	typeRNum *pmax_ = grampc->pmax;

	if (grampc->ScaleProblem == INT_ON) {
		pmin_ = grampc->rwsScale;
		pmax_ = pmin_ + grampc->Np;
		scale_parameters(pmin_, grampc->pmin, grampc);
		scale_parameters(pmax_, grampc->pmax, grampc);
	}

	for (i = 0; i < grampc->Np; i++) {
		// lower bound
		if (p[i] < pmin_[i]) {
			p[i] = pmin_[i];
		}
		// upper bound
		else if (p[i] > pmax_[i]) {
			p[i] = pmax_[i];
		}
	}
}
*/
void evaluate_gradT(typeGRAMPC *grampc)
{
	typeInt i, k;

	ctypeRNum *t = grampc->t;
	ctypeRNum *x_ = grampc->x + (grampc->Nhor - 1) * grampc->Nx;
	ctypeRNum *adj_ = grampc->adj + (grampc->Nhor - 1) * grampc->Nx;
	ctypeRNum *u_ = grampc->u + (grampc->Nhor - 1) * grampc->Nu;
	typeRNum *p_ = grampc->p;

	ctypeRNum *mult = grampc->mult + (grampc->Nhor - 1) * grampc->Nc;
	ctypeRNum *pen = grampc->rws_pen + (grampc->Nhor - 1) * grampc->Nc;
	ctypeRNum *cfct = grampc->rws_cfct + (grampc->Nhor - 1) * grampc->Nc;

	typeRNum *s1 = grampc->rwsGeneral;
	typeRNum l = 0.0;
	typeRNum f = 0.0;
	typeRNum c = 0.0;
	typeRNum dVdt = 0.0;

	/* Unscaling */
	if (grampc->ScaleProblem == INT_ON) {
		ASSIGN_X(x_, x_, grampc);
		ASSIGN_ADJ(adj_, adj_, grampc);
		ASSIGN_U(u_, u_, grampc);
		ASSIGN_P(p_, p_, grampc);
	}

	/* H(x(T),u(T),p,adj(T),T) = */
	/* l(x,u,p,t)  */
	if (grampc->IntegralCost == INT_ON) {
		lfct(&l, t[grampc->Nhor - 1], x_, u_, p_, grampc->xdes, grampc->udes, grampc->userparam);
/*		if (grampc->ScaleProblem == INT_ON) {
			scale_cost(&l, grampc->JScale, 1);
		}   */
	}

	/* + adj' * f(x,u,p,t)   */
//	ffct(s1, t[grampc->Nhor - 1], x_, u_, p_, grampc->userparam);
ffct(s1, t[grampc->Nhor - 1], x_, u_, grampc->userparam);
	MatMult(&f, adj_, s1, 1, grampc->Nx, 1);

	/* + c(x,u,p,t)' * (mult + c(x,u,p,t)' * diag(pen) / 2); */
	if (grampc->EqualityConstraints == INT_ON) {
		for (i = 0; i < grampc->Ng; i++) {
			k = i;
			c = c + cfct[k] * (mult[k] + cfct[k] * pen[k] / 2);
		}
	}
	if (grampc->InequalityConstraints == INT_ON) {
		for (i = 0; i < grampc->Nh; i++) {
			k = i + grampc->Ng;
			c = c + cfct[k] * (mult[k] + cfct[k] * pen[k] / 2);
		}
	}

	/* dV(x,p,T)/dT */
	if (grampc->TerminalCost == INT_ON) {
		dVdT(&dVdt, t[grampc->Nhor - 1], x_, p_, grampc->xdes, grampc->userparam);
/*		if (grampc->ScaleProblem == INT_ON) {
			scale_cost(&dVdt, grampc->JScale, 1);
		}   */
	}

	/* + dc(x,p,T)/dT' * mult(T) */
	grampc->gradT = l + f + c + dVdt + grampc->dcdt;

	/* Scaling */
	if (grampc->ScaleProblem == INT_ON) {
		grampc->gradT = grampc->gradT * grampc->TScale;
	}
}

void timeproj(typeRNum *T, typeGRAMPC *grampc)
{
	typeRNum Tmin_ = grampc->Tmin;
	typeRNum Tmax_ = grampc->Tmax;

	if (grampc->ScaleProblem == INT_ON) {
		scale_time(&Tmin_, Tmin_, grampc);
		scale_time(&Tmax_, Tmax_, grampc);
	}

	/* lower bound */
	if (T[0] < Tmin_) {
		T[0] = Tmin_;
	}
	/* upper bound */
	else if (T[0] > Tmax_) {
		T[0] = Tmax_;
	}
}
/*
void linesearch_adaptive(typeRNum *alpha, ctypeInt igrad, typeGRAMPC *grampc)
{
	typeInt ils, j;
	ctypeInt NLS_adapt = 2 * (NALS + 1);
	typeRNum Jls[2] = { 0 , 0 };

	typeRNum *t = grampc->t;
	typeRNum *tls = grampc->tls;
	typeRNum *u = grampc->u;
	typeRNum *uls = grampc->uls;
	ctypeRNum *gradu = grampc->gradu;
	typeRNum *p = grampc->p;
	typeRNum *pls = grampc->pls;
//	ctypeRNum *gradp = grampc->gradp;
	ctypeRNum *T = &grampc->T;
	typeRNum Tls;
	ctypeRNum *gradT = &grampc->gradT;
	typeRNum *lsAdapt = grampc->lsAdapt + grampc->MaxGradIter*NLS_adapt;

	// Adapt interval
	if (ABS(lsAdapt[NALS + 1] - lsAdapt[NALS + 3]) > grampc->LineSearchAdaptAbsTol) {
		if (lsAdapt[NALS] >= lsAdapt[0] + (1 - grampc->LineSearchIntervalTol)*(lsAdapt[NALS - 1] - lsAdapt[0])) {
			if (lsAdapt[NALS - 1] <= grampc->LineSearchMax) {
				for (ils = 0; ils < NALS; ils++) {
					lsAdapt[ils] = lsAdapt[ils] * grampc->LineSearchAdaptFactor;
				}
			}
			else {
				//grampc->status |= STATUS_LINESEARCH_MAX;
			}
		}
		else if (lsAdapt[NALS] <= lsAdapt[0] + grampc->LineSearchIntervalTol*(lsAdapt[NALS - 1] - lsAdapt[0])) {
			if (lsAdapt[0] >= grampc->LineSearchMin) {
				for (ils = 0; ils < NALS; ils++) {
					lsAdapt[ils] = lsAdapt[ils] / grampc->LineSearchAdaptFactor;
				}
			}
			else {
				//grampc->status |= STATUS_LINESEARCH_MIN;
			}
		}
	}

	// Evaluate cost for NLS step sizes
	for (ils = 0; ils < NALS; ils++) {

		if (grampc->OptimControl == INT_ON) {
			for (j = 0; j < grampc->Nhor*grampc->Nu; j++) {
				uls[j] = u[j] - lsAdapt[ils] * gradu[j];
			}
			inputproj(uls, grampc);
		}
		else {
			uls = u;
		}

 		if (grampc->OptimParam == INT_ON) {
			for (j = 0; j < grampc->Np; j++) {
				pls[j] = p[j] - grampc->OptimParamLineSearchFactor * lsAdapt[ils] * gradp[j];
			}
			paramproj(pls, grampc);
		}
		 else {
			pls = p;
	 	}

		if (grampc->OptimTime == INT_ON) {
			Tls = T[0] - grampc->OptimTimeLineSearchFactor * lsAdapt[ils] * gradT[0];
			timeproj(&Tls, grampc);
			discretize_time(tls, Tls, &grampc);
		}
		else {
			tls = t;
		}

		evaluate_sys(tls, uls, pls, grampc);
		evaluate_constraints(tls, uls, pls, 0, 0, grampc);
		evaluate_cost(Jls, tls, uls, pls, grampc);
		lsAdapt[ils + NALS + 1] = Jls[1];
	}

	// Determine optimal step size by curve fitting
	lsearch_fit(lsAdapt + NALS, lsAdapt + 2 * NALS + 1, lsAdapt, lsAdapt + NALS + 1);
	*alpha = lsAdapt[NALS];

	// Save history
	for (ils = 0; ils < 2 * (NALS + 1); ils++) {
		grampc->lsAdapt[ils + igrad * NLS_adapt] = lsAdapt[ils];
	}

}
*/
void linesearch_explicit(typeRNum *alpha, typeGRAMPC *grampc)
{
	typeInt i, j;
	typeRNum *lsExplicit = grampc->lsExplicit;
	typeRNum graduMax, graduAbs;
//	typeRNum NomDenum[2];

	if (grampc->OptimControl == INT_ON) {
		update_lsExplicit(lsExplicit, grampc->u, grampc->uprev, grampc->gradu, grampc->graduprev, grampc->Nhor*grampc->Nu, grampc);
	}
	//OptimParam=off
/*	if (grampc->OptimParam == INT_ON) {
		update_lsExplicit(NomDenum, grampc->p, grampc->pprev, grampc->gradp, grampc->gradpprev, grampc->Np, grampc);
		lsExplicit[0] += NomDenum[0] * grampc->OptimParamLineSearchFactor;
		lsExplicit[1] += NomDenum[1] * grampc->OptimParamLineSearchFactor * grampc->OptimParamLineSearchFactor;
	}*/

	//OptimTime=off
/*	if (grampc->OptimTime == INT_ON) {
		update_lsExplicit(NomDenum, &grampc->T, &grampc->Tprev, &grampc->gradT, &grampc->gradTprev, 1, grampc);
		lsExplicit[0] += NomDenum[0] * grampc->OptimTimeLineSearchFactor;
		lsExplicit[1] += NomDenum[1] * grampc->OptimTimeLineSearchFactor * grampc->OptimTimeLineSearchFactor;
	}*/

	if (lsExplicit[0] > 0 && lsExplicit[1] > 0) {
		lsExplicit[2] = lsExplicit[0] / lsExplicit[1];
	}
	else {
		/* AutoFallback Method if the option is set OptimControl is on and limits for every input are set */
		if (grampc->LineSearchExpAutoFallback == INT_ON && grampc->OptimControl == INT_ON && lsExplicit[3] == 1) {
			/* Limit the stepsize to 10% of the maximum value */
			lsExplicit[2] = grampc->LineSearchMax / 10;

			/* Set the stepsize, that the control update is not more than 1.0% of the control range */
			for (j = 0; j < grampc->Nu; j++) {
				graduMax = 0;
				for (i = 0; i < grampc->Nhor; i++) {
					graduAbs = ABS(grampc->gradu[j + i * grampc->Nu]);
					if (graduAbs > graduMax) {
						graduMax = graduAbs;
					}
				}
				if (grampc->ScaleProblem == INT_ON) {
					graduMax = graduMax / grampc->uScale[j];
				}
				lsExplicit[2] = MIN(lsExplicit[2], (grampc->umax[j] - grampc->umin[j]) / (graduMax * 100));
			}
		}
		else {
			lsExplicit[2] = grampc->LineSearchInit;
			//grampc->status |= STATUS_LINESEARCH_INIT;
		}
	}

	if (lsExplicit[2] > grampc->LineSearchMax) {
		lsExplicit[2] = grampc->LineSearchMax;
		//grampc->status |= STATUS_LINESEARCH_MAX;
	}
	else if (lsExplicit[2] < grampc->LineSearchMin) {
		lsExplicit[2] = grampc->LineSearchMin;
		//grampc->status |= STATUS_LINESEARCH_MIN;
	}

	*alpha = lsExplicit[2];
}

void update_lsExplicit(typeRNum* NomDenum, ctypeRNum* a, ctypeRNum* aprev, ctypeRNum* dHda, ctypeRNum* dHdaprev, ctypeInt length, typeGRAMPC *grampc) {
	typeRNum diffa, diffdHda;
	typeInt j;

	NomDenum[0] = 0;
	NomDenum[1] = 0;

	if (grampc->LineSearchType == INT_EXPLS1) {
		/* Formula #1 */
		for (j = 0; j < length; j++) {
			diffa = a[j] - aprev[j];
			diffdHda = dHda[j] - dHdaprev[j];
			NomDenum[0] = NomDenum[0] + diffa * diffa;
			NomDenum[1] = NomDenum[1] + diffa * diffdHda;
		}
	}
	else {
		/* Formula #2 */
		for (j = 0; j < length; j++) {
			diffa = a[j] - aprev[j];
			diffdHda = dHda[j] - dHdaprev[j];
			NomDenum[0] = NomDenum[0] + diffa * diffdHda;
			NomDenum[1] = NomDenum[1] + diffdHda * diffdHda;
		}
	}
}

void evaluate_cost(typeRNum *s, ctypeRNum *t, ctypeRNum *u, ctypeRNum *p, typeGRAMPC *grampc)
{
	//typeInVfctPtr pIntCost;
	typeRNum Jint[2] = { 0, 0 };
	typeRNum Jterm[2] = { 0, 0 };
	ctypeRNum *p_ = p;

	/* integrator for integral cost function */
/*	if (grampc->IntegratorCost == INT_TRAPZ) {
		pIntCost = &trapezodial;
	}
	else {
		pIntCost = &simpson;
	}
*/
	/* Unscaling */
	if (grampc->ScaleProblem == INT_ON) {
		ASSIGN_P(p_, p, grampc);
	}

	/* integrate cost */
	trapezodial(Jint, t, grampc->x, u, p_, grampc);

	/* terminal cost */
	WtermCost(Jterm, t[grampc->Nhor - 1],
		grampc->x + (grampc->Nhor - 1) * grampc->Nx, p_,
		grampc->mult + (grampc->Nhor - 1) * grampc->Nc,
		grampc->rws_pen + (grampc->Nhor - 1) * grampc->Nc,
		grampc->rws_cfct + (grampc->Nhor - 1) * grampc->Nc, grampc);

	s[0] = Jint[0] + Jterm[0];
	s[1] = Jint[1] + Jterm[1];
}

void WintCost(typeRNum *s, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p_,
	ctypeRNum *mult, ctypeRNum *pen, ctypeRNum *cfct, typeGRAMPC *grampc)
{
	typeInt i, k;
	ctypeRNum *x_ = x;
	ctypeRNum *u_ = u;

	/* Unscaling */
	if (grampc->ScaleProblem == INT_ON) {
		ASSIGN_X(x_, x, grampc);
		ASSIGN_U(u_, u, grampc);
	}

	s[0] = 0;
	s[1] = 0;

	/* Integral Cost */
	if (grampc->IntegralCost == INT_ON) {
		lfct(s, t, x_, u_, p_, grampc->xdes, grampc->udes, grampc->userparam);
/*		if (grampc->ScaleProblem == INT_ON) {
			scale_cost(s, grampc->JScale, 1);
		}*/
	}

	/* Equality Constraints */
	if (grampc->EqualityConstraints == INT_ON) {
		for (i = 0; i < grampc->Ng; i++) {
			k = i;
			s[1] = cfct[k] * (mult[k] + cfct[k] * pen[k] / 2);
		}
	}

	/* Inequality Constraints */
	if (grampc->InequalityConstraints == INT_ON) {
		for (i = 0; i < grampc->Nh; i++) {
			k = grampc->Ng + i;
			s[1] = s[1] + cfct[k] * (mult[k] + cfct[k] * pen[k] / 2);
		}
	}
	s[1] = s[1] + s[0];
}

void WtermCost(typeRNum *s, ctypeRNum t, ctypeRNum *x, ctypeRNum *p_,
	ctypeRNum *mult, ctypeRNum *pen, ctypeRNum *cfct, typeGRAMPC *grampc)
{
	typeInt i, k;
	ctypeRNum *x_ = x;

	/* Unscaling */
	if (grampc->ScaleProblem == INT_ON) {
		ASSIGN_X(x_, x, grampc);
	}

	s[0] = 0;
	s[1] = 0;

	/* Terminal Cost */
	if (grampc->TerminalCost == INT_ON) {
		Vfct(s, t, x_, p_, grampc->xdes, grampc->userparam);
/*		if (grampc->ScaleProblem == INT_ON) {
			scale_cost(s, grampc->JScale, 1);
		}*/
	}

	/* Terminal Equality Constraints */
	if (grampc->TerminalEqualityConstraints == INT_ON) {
		for (i = 0; i < grampc->NgT; i++) {
			k = grampc->Ng + grampc->Nh + i;
			s[1] = cfct[k] * (mult[k] + cfct[k] * pen[k] / 2);
		}
	}

	/* Terminal Inequality Constraints */
	if (grampc->TerminalInequalityConstraints == INT_ON) {
		for (i = 0; i < grampc->NhT; i++) {
			k = grampc->Ng + grampc->Nh + grampc->NgT + i;
			s[1] = s[1] + cfct[k] * (mult[k] + cfct[k] * pen[k] / 2);
		}
	}
	s[1] = s[1] + s[0];
}

typeBoolean convergence_test_gradient(typeRNum rel_tol, typeGRAMPC *grampc)
{
	typeRNum u_diff_norm = 0;
	typeRNum u_norm = 0;
	typeRNum u_rel_norm = 0;
//	typeRNum p_diff_norm = 0;
//	typeRNum p_norm = 0;
	typeRNum p_rel_norm = 0;
	typeRNum T_diff_norm = 0;
	typeRNum T_norm = 0;
	typeRNum T_rel_norm = 0;

	if (grampc->OptimControl == INT_ON) {
		/* ||u - uprev|| */
		MatDiffNorm(&u_diff_norm, grampc->u, grampc->uprev, grampc->Nhor, grampc->Nu);
		/* ||u|| */
		MatNorm(&u_norm, grampc->u, grampc->Nhor, grampc->Nu);
		u_rel_norm = u_norm > 0 ? u_diff_norm / u_norm : 0;
	}
/*	if (grampc->OptimParam == INT_ON) {
		// ||p - pprev||
		MatDiffNorm(&p_diff_norm, grampc->p, grampc->pprev, 1, grampc->Np);
		// ||p||
		MatNorm(&p_norm, grampc->p, 1, grampc->Np);
		p_rel_norm = p_norm > 0 ? p_diff_norm / p_norm : 0;
	}*/
	if (grampc->OptimTime == INT_ON) {
		/* ||T - Tprev|| */
		T_diff_norm = (grampc->T - grampc->Tprev) * (grampc->T - grampc->Tprev);
		/* ||T|| */
		T_norm = grampc->T * grampc->T;
		T_rel_norm = T_norm > 0 ? SQRT(T_diff_norm / T_norm) : 0;
	}
	return MAX(MAX(u_rel_norm, p_rel_norm), T_rel_norm) < rel_tol;
}

typeBoolean convergence_test_constraints(ctypeRNum *abs_tol, typeGRAMPC *grampc)
{
	typeInt i, j, k;
	ctypeRNum *cfct = grampc->rws_cfct;

	if (grampc->EqualityConstraints == INT_ON) {
		for (j = 0; j < grampc->Ng; j++) {
			k = j;
			for (i = 1; i < grampc->Nhor; i++) {
				if (ABS(cfct[i * grampc->Nc + k]) > abs_tol[k]) {
					return 0;
				}
			}
		}
	}
	if (grampc->InequalityConstraints == INT_ON) {
		for (j = 0; j < grampc->Nh; j++) {
			k = grampc->Ng + j;
			for (i = 1; i < grampc->Nhor; i++) {
				if (cfct[i * grampc->Nc + k] > abs_tol[k]) {
					return 0;
				}
			}
		}
	}

	i = grampc->Nhor - 1;
	if (grampc->TerminalEqualityConstraints == INT_ON) {
		for (j = 0; j < grampc->NgT; j++) {
			k = grampc->Ng + grampc->Nh + j;
			if (ABS(cfct[i * grampc->Nc + k]) > abs_tol[k]) {
				return 0;
			}
		}
	}
	if (grampc->TerminalInequalityConstraints == INT_ON) {
		for (j = 0; j < grampc->NhT; j++) {
			k = grampc->Ng + grampc->Nh + grampc->NgT + j;
			if (cfct[i * grampc->Nc + k] > abs_tol[k]) {
				return 0;
			}
		}
	}
	return 1;
}

void shiftTrajectory(typeRNum *trajectory, ctypeInt Nhor, ctypeInt Nrows, ctypeInt Nshiftrows, ctypeRNum dt, ctypeRNum *t)
{
	/* function assumes uniform grid size over prediction horizon */
	typeInt i, j;

	/* compute how far the sampling points must be shifted */
	ctypeInt shift = (typeInt)(dt / (t[1] - t[0]));
	ctypeInt shiftvar = shift * Nrows;
	ctypeRNum interpfact = (dt / (t[1] - t[0])) - shift;

/*	if (shift >= Nhor) {
		grampc_error("Horizon too short for the current sampling time.");
	}*/

	/* Interpolation between the grid points */
	for (i = 0; i < Nhor - 1 - shift; i++) {
		for (j = 0; j < Nshiftrows; j++) {
			trajectory[j] = trajectory[shiftvar + j] + (trajectory[shiftvar + j + Nrows] - trajectory[shiftvar + j]) * interpfact;
		}
		trajectory += Nrows;
	}

	/* next element: extrapolation */
 /*
	if (shift == Nhor - 1) {
		for (j = 0; j < Nrows; j++) {
			trajectory[j] = trajectory[shiftvar + j] + (trajectory[shiftvar + j] - trajectory[shiftvar + j - Nrows]) * interpfact;
		}
	}
	else {
		for (j = 0; j < Nrows; j++) {
			trajectory[j] = trajectory[shiftvar + j] + (trajectory[shiftvar + j] - trajectory[j - Nrows]) / (1 / interpfact - 1);
		}
	}
	trajectory += Nrows;
	i++;
	*/

	/* if there are elements left hold last value */
	for (; i < Nhor; i++) {
		for (j = 0; j < Nshiftrows; j++) {
			trajectory[j] = trajectory[j - Nrows];
		}
		trajectory += Nrows;
	}
}

void shortenTrajectory(typeRNum *trajectory, ctypeInt Nhor, ctypeInt Nrows, ctypeInt Nshortenrows, ctypeRNum dt, ctypeRNum *t)
{
	/* function assumes uniform grid size over prediction horizon */
	typeInt i, j, ind;
	typeRNum takt;
	typeRNum dtgn = (t[Nhor - 1] - dt) / (Nhor - 1); /* compute new grid step size */

	/* Interpolation between the grid points from the beginning to the end */
	for (i = 0; i < Nhor - 1; i++) {
		takt = dt + dtgn * i;
		ind = (typeInt)(takt / (t[1] - t[0]));

		for (j = 0; j < Nshortenrows; j++) {
			trajectory[j + i * Nrows] = trajectory[ind*Nrows + j] + (trajectory[(ind + 1)* Nrows + j] - trajectory[ind*Nrows + j]) *(takt - t[ind]) / (t[1] - t[0]);
		}
	}
}

void intsysEuler_wsys(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum * t, ctypeRNum * x,
	ctypeRNum *u, ctypeRNum *p_, typeGRAMPC *grampc)
{
	typeInt i, j;
	typeRNum h;
	ctypeRNum *dcdx = grampc->dcdx + grampc->Nx * (grampc->Nhor - 1);
	typeRNum *s = grampc->rwsGeneral + LWadjsys; /* size Nx */

	for (j = 0; j < Nint - 1; j++) {
		if (j > 0) {
			t += pInt;
			x += pInt * grampc->Nx;
			u += pInt * grampc->Nu;
			y += pInt * grampc->Nx;

			dcdx += (-1)*grampc->Nx; /* only used for integrating the adjoint system */
		}

	//	(*pfct)(s, y, t, x, u, p_, dcdx, grampc);
	Wsys(s, y, t, x, u, p_, dcdx, grampc);

		h = t[pInt] - t[0];
		for (i = 0; i < grampc->Nx; i++) {
			y[i + pInt * grampc->Nx] = y[i] + h * s[i];
		}
	}
}

void intsysEuler_wadjsys(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum * t, ctypeRNum * x,
	ctypeRNum *u, ctypeRNum *p_, typeGRAMPC *grampc)
{
	typeInt i, j;
	typeRNum h;
	ctypeRNum *dcdx = grampc->dcdx + grampc->Nx * (grampc->Nhor - 1);
	typeRNum *s = grampc->rwsGeneral + LWadjsys; /* size Nx */

	for (j = 0; j < Nint - 1; j++) {
		if (j > 0) {
			t += pInt;
			x += pInt * grampc->Nx;
			u += pInt * grampc->Nu;
			y += pInt * grampc->Nx;

			dcdx += (-1)*grampc->Nx; /* only used for integrating the adjoint system */
		}

	//	(*pfct)(s, y, t, x, u, p_, dcdx, grampc);
	Wadjsys(s, y, t, x, u, p_, dcdx, grampc);

		h = t[pInt] - t[0];
		for (i = 0; i < grampc->Nx; i++) {
			y[i + pInt * grampc->Nx] = y[i] + h * s[i];
		}
	}
}

void trapezodial(typeRNum *s, ctypeRNum *t, ctypeRNum *x, ctypeRNum *u,
	ctypeRNum *p_, typeGRAMPC *grampc)
{
	typeInt i;
	typeRNum h;

	ctypeRNum *mult = grampc->mult;
	ctypeRNum *pen = grampc->rws_pen;
	ctypeRNum *cfct = grampc->rws_cfct;

	typeRNum *s1 = grampc->rwsGeneral; /* size: 2*/

	s[0] = 0;
	s[1] = 0;

	/* Integration */
	for (i = 0; i < grampc->Nhor; i++) {

		WintCost(s1, t[i], x + i * grampc->Nx, u + i * grampc->Nu, p_,
			mult + i * grampc->Nc, pen + i * grampc->Nc, cfct + i * grampc->Nc, grampc);

		if (i == 0) {
			h = (t[i + 1] - t[i]) / 2;
		}
		else if (i <= grampc->Nhor - 2) {
			h = (t[i + 1] - t[i - 1]) / 2;
		}
		else {
			h = (t[i] - t[i - 1]) / 2;
		}

		s[0] = s[0] + h * s1[0];
		s[1] = s[1] + h * s1[1];
	}
}

#ifdef __cplusplus
}
#endif // __cplusplus
