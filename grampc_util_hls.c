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



#include "grampc_util.h"
#include "grampc_run.h"
/*
void createNumMatrix(typeRNum **cs, const size_t size) {
	if (size == 0) {
		*cs = NULL;
	}
	else {
		*cs = (typeRNum *)calloc(size, sizeof(typeRNum));
		if (*cs == NULL) {
			grampc_error(PARAM_ALLOC_FAILED);
		}
	}
}

void createIntMatrix(typeInt **cs, const size_t size) {
	if (size == 0) {
		*cs = NULL;
	}
	else {
		*cs = (typeInt *)calloc(size, sizeof(typeInt));
		if (*cs == NULL) {
			grampc_error(PARAM_ALLOC_FAILED);
		}
	}
}
*/
/*
void resizeNumMatrix(typeRNum **cs, const size_t size) {
	free(*cs);
	if (size == 0) {
		*cs = NULL;
	}
	else {
		*cs = (typeRNum *)calloc(size, sizeof(typeRNum));
		if (*cs == NULL) {
			grampc_error(RWS_ELEMENT_ALLOC_FAILED);
		}
	}
}

void resizeIntMatrix(typeInt **cs, const size_t size) {
	free(*cs);
	if (size == 0) {
		*cs = NULL;
	}
	else {
		*cs = (typeInt *)calloc(size, sizeof(typeInt));
		if (*cs == NULL) {
			grampc_error(SOL_ALLOC_FAILED);
		}
	}
}
*/
void init_rws_time(typeGRAMPC **grampc) {
/*	if (grampc->ScaleProblem == INT_ON) {
		scale_time(&grampc->T, grampc->Thor, grampc);
	}
	else {*/
		(*grampc)->T = (*grampc)->Thor;
	//}
	(*grampc)->Tprev = (*grampc)->T;
	//discretize_time((*grampc)->t, (*grampc)->T, &(*grampc));
	char i;
	char N = (*grampc)->Nhor;
	typeRNum m, c;
	typeRNum T_ = (*grampc)->T;
	typeRNum *tvec = (*grampc)->t;

	if ((*grampc)->TimeDiscretization == INT_UNIFORM || T_ <= (N - 1) * (*grampc)->dt) {
            m = T_ / (N - 1);
            for (i = 0; i < (*grampc)->Nhor; i++) {
                tvec[i] = m * i;
                }
	}
	/* Nonuniform discretization of interval [0, T_] */
	else if ((*grampc)->TimeDiscretization == INT_NONUNIFORM) {
		m = (T_ / (N - 1) - (*grampc)->dt) / (N - 2);
		c = (*grampc)->dt - m;
		for (i = 0; i < (*grampc)->Nhor; i++) {
			tvec[i] = m * i * i + c * i;
		}
	}
}

void init_rws_controls(typeGRAMPC **grampc) {
	char i, k;
	for (i = 0; i < (*grampc)->Nhor; i++) {
		k = i * (*grampc)->Nu;
/*		if (grampc->ScaleProblem == INT_ON) {
			scale_controls(grampc->u + k, grampc->u0, grampc);
		}
		else {*/
			MatCopy((*grampc)->u + k, (*grampc)->u0, 1, (*grampc)->Nu);
            MatCopy((*grampc)->uprev, (*grampc)->u, (*grampc)->Nhor, (*grampc)->Nu);
   //     }
	}
}

void resize_rwsGeneral(typeGRAMPC **grampc) {
	typeInt LInt = LWadjsys;
	typeInt LIntCost = 0;
	typeInt LConst = 0;
/*
	switch (grampc->Integrator) {
	case INT_EULER:    LInt += Leuler; break;
	case INT_MODEULER: LInt += Lmodeuler; break;
	case INT_HEUN:     LInt += Lheun; break;
	case INT_RODAS:		 LInt += Lrodas; break;
	case INT_RUKU45:   LInt += Lruku45; break;
	}
*/
	LInt += Leuler;
/*
	switch (grampc->IntegratorCost) {
	case INT_TRAPZ:   LIntCost = LIntCostTrapezoidal; break;
	case INT_SIMPSON: LIntCost = LIntCostSimpson; break;
	}*/

	LIntCost = LIntCostTrapezoidal;

	if ((*grampc)->Nc > 0) {
		LConst = LevaluateConstraints;
	}

	(*grampc)->lrwsGeneral = MAX(MAX(MAX(MAX(MAX(Lgradp, Lgradu), LgradT), LIntCost), LInt), LConst);
	//resizeNumMatrix(&grampc->rwsGeneral, grampc->lrwsGeneral);
}
/*
void resize_rwsRodas(typeGRAMPC *grampc) {
	if (grampc->Integrator == INT_RODAS) {
		resizeNumMatrix(&grampc->rparRodas, grampc->Nx*grampc->Nhor);
		resizeIntMatrix(&grampc->iparRodas, 20);
		resizeNumMatrix(&grampc->workRodas, grampc->lworkRodas);
		resizeIntMatrix(&grampc->iworkRodas, grampc->liworkRodas);
	}
	else {
		resizeNumMatrix(&grampc->rparRodas, 0);
		resizeIntMatrix(&grampc->iparRodas, 0);
		resizeNumMatrix(&grampc->workRodas, 0);
		resizeIntMatrix(&grampc->iworkRodas, 0);
	}
}


void setLWorkRodas(typeGRAMPC *grampc)
{
	typeInt LJAC, LMAS, LE1;												// length of rodas work space
																			// ctypeInt IFCN = grampc->FlagsRodas[0]; 		  0 --> right hand side independent of time t
																			// ctypeInt IDFX = grampc->FlagsRodas[1];		  0 --> DF/DX is numerically computed
																			// typeInt IJAC = grampc->FlagsRodas[2];		  1(0) -> analytical (numerical) jacobian (partial derivatives of right hand side w.r.t. state)
	ctypeInt MLJAC = grampc->FlagsRodas[4];		// no. of lower diagonals of jacobian
	ctypeInt MUJAC = grampc->FlagsRodas[5];		// no. of upper diagonals of jacobian

	typeInt IMAS = grampc->FlagsRodas[3];		// 1 --> mass matrix is supplied
	ctypeInt MLMAS = grampc->FlagsRodas[6];		// no. of lower diagonals of mass matrix
	ctypeInt MUMAS = grampc->FlagsRodas[7];		// no. of upper diagonals of mass matrix

	if (MLJAC < grampc->Nx)
	{
		LJAC = MLJAC + MUJAC + 1;
		LE1 = 2 * MLJAC + MUJAC + 1;
	}
	else
	{
		LJAC = grampc->Nx;
		LE1 = grampc->Nx;
	}

	if (IMAS == 0)
	{
		LMAS = 0;
	}
	else
	{
		if (MLMAS == grampc->Nx)
		{
			LMAS = grampc->Nx;
		}
		else
		{
			LMAS = MLMAS + MUMAS + 1;
		}
	}
	grampc->lworkRodas = grampc->Nx*(LJAC + LMAS + LE1 + 14) + 20;
//	resizeNumMatrix(&grampc->workRodas, grampc->lworkRodas);
}
*/
/*
typeInt CastDvec2Intvec(typeInt** Intvec, const double* Numvec, const size_t size) {
	unsigned typeInt i;
	*Intvec = (typeInt*)malloc(size * sizeof(typeInt));
	if (*Intvec != NULL) {
		for (i = 0; i < size; i++) {
			(*Intvec)[i] = (typeInt)Numvec[i];
		}
		return 1;
	}
	return -1;
}

typeInt CastDvec2Numvec(typeRNum** Realvec, const double* Numvec, const size_t size) {
	unsigned typeInt i;
	*Realvec = (typeRNum*)malloc(size * sizeof(typeRNum));
	if (*Realvec != NULL) {
		for (i = 0; i < size; i++) {
			(*Realvec)[i] = (typeRNum)Numvec[i];
		}
		return 1;
	}
	return -1;
}
*/

void discretize_time(typeRNum *tvec, typeRNum T, typeGRAMPC **grampc)
{
	char i;
	char N = (*grampc)->Nhor;
	typeRNum m, c;
	typeRNum T_ = T;

/*	if (grampc->ScaleProblem == INT_ON) {
		unscale_time(&T_, T, grampc);
	}*/

	/* Uniform discretization of interval [0, T_] */
	if ((*grampc)->TimeDiscretization == INT_UNIFORM || T_ <= (N - 1) * (*grampc)->dt) {
		m = T_ / (N - 1);
		for (i = 0; i < (*grampc)->Nhor; i++) {
			tvec[i] = m * i;
		}
	}
	/* Nonuniform discretization of interval [0, T_] */
	else if ((*grampc)->TimeDiscretization == INT_NONUNIFORM) {
		m = (T_ / (N - 1) - (*grampc)->dt) / (N - 2);
		c = (*grampc)->dt - m;
		for (i = 0; i < (*grampc)->Nhor; i++) {
			tvec[i] = m * i * i + c * i;
		}
	}
}
/*
void check_ControlLimits(typeGRAMPC* grampc) {
	typeInt i;
	// if any control limit is INF or -INF disable auto fallback strategy of the explicit linesearch */
/*	if (grampc->LineSearchType == INT_EXPLS1 || grampc->LineSearchType == INT_EXPLS2) {
		grampc->lsExplicit[3] = 1;

		for (i = 0; i < grampc->Nu; i++) {
			if (grampc->umax[i] >= INF || grampc->umin[i] <= -INF) {
				grampc->lsExplicit[3] = 0;
			}
		}
	}
}
*/
char* IntegratorInt2Str(ctypeInt INT_Integrator) {
	switch (INT_Integrator) {
	case INT_EULER:    return "euler";
	case INT_MODEULER: return "modeuler";
	case INT_HEUN:     return "heun";
	case INT_RODAS:    return "rodas";
	case INT_RUKU45:   return "ruku45";
	}
	return "";
}

char* LineSearchTypeInt2Str(ctypeInt INT_LineSearchType) {
	switch (INT_LineSearchType) {
	case INT_ADAPTIVELS: return "adaptive";
	case INT_EXPLS1:     return "explicit1";
	case INT_EXPLS2:     return "explicit2";
	}
	return "";
}


void lsearch_fit(typeRNum *kfit, typeRNum *Jfit, ctypeRNum *k, ctypeRNum *J)
{
	ctypeRNum aux = (2 * J[1] - J[0] - J[2]);
	ctypeRNum aux2 = (J[2] - 4 * J[1]);

	/* exists a minimum? (positive curvature) */
	if (aux <= -aEPS) {
		kfit[0] = k[1] + (k[1] - k[0]) / 2 * (J[2] - J[0]) / aux;
		Jfit[0] = (J[0] * J[0] + aux2 * aux2 - 2 * J[0] * (4 * J[1] + J[2])) / (8 * aux);
	}
	/* smallest J */
	if (aux > -aEPS  || kfit[0] < k[0] || kfit[0] > k[2]) {
		if (J[0] <= J[1] && J[0] <= J[2]) {
			kfit[0] = k[0];
			Jfit[0] = J[0];
		}
		else if (J[2] <= J[0]  && J[2] <= J[1]) {
			kfit[0] = k[2];
			Jfit[0] = J[2];
		}
		else {
			kfit[0] = k[1];
			Jfit[0] = J[1];
		}
	}
}

void interplin(typeRNum *varint, ctypeRNum *tvec, ctypeRNum *varvec, ctypeRNum tint,
	ctypeInt Nvar, ctypeInt Nvec, ctypeInt searchdir)
{
	/* option to determine position ioff in time vector such that tvec[ioff] <= tint <= tvec[ioff+1]
	 * searchdir= 1: search in forward direction starting at ioff=0									 *
	 * searchdir=-1: search in backward direction starting at ioff=Nvec-1
	 */
	typeInt i;
	typeInt ioff;
	typeRNum dtratio;
	ctypeRNum *var0, *var1;

	if (tint <= tvec[0]) {
		for (i = 0; i < Nvar; i++) {
			varint[i] = varvec[i];
		}
	}
	else if (tint >= tvec[Nvec - 1]) {
		var0 = varvec + (Nvec - 1)*Nvar;
		for (i = 0; i < Nvar; i++)
			varint[i] = var0[i];
	}
	else {
		/* ---- HLS PATCH ----------------------------------------------------
		 * The original used two UNBOUNDED while loops here:
		 *     if (searchdir==1){ioff=0;        while(tvec[ioff]<tint) ioff+=1; ioff-=1;}
		 *     else             {ioff=Nvec-2;   while(tvec[ioff]>tint) ioff-=1;}
		 * They terminate in practice (bounded by Nvec = Nhor = 19) but HLS
		 * cannot prove it, so it cannot schedule them. Rewritten as bounded
		 * for-loops with an explicit trip count + break. Semantics identical.
		 * Note the two branches ABOVE already handled tint<=tvec[0] and
		 * tint>=tvec[Nvec-1], so here tvec[0] < tint < tvec[Nvec-1] and a valid
		 * ioff always exists; the clamps just make that provable.
		 * ------------------------------------------------------------------ */
		if (searchdir == 1) {
			ioff = 0;
			for (i = 0; i < Nvec; i++) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=19
				if (tvec[i] >= tint) { ioff = i; break; }
				ioff = i + 1;
			}
			ioff -= 1;
			if (ioff < 0)          ioff = 0;
			if (ioff > Nvec - 2)   ioff = Nvec - 2;
		}
		else {
			ioff = Nvec - 2;
			for (i = Nvec - 2; i >= 0; i--) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=19
				if (tvec[i] <= tint) { ioff = i; break; }
				ioff = i - 1;
			}
			if (ioff < 0)          ioff = 0;
			if (ioff > Nvec - 2)   ioff = Nvec - 2;
		}
		dtratio = (tint - tvec[ioff]) / (tvec[ioff + 1] - tvec[ioff]);
		var0 = varvec + ioff * Nvar;
		var1 = var0 + Nvar;
		for (i = 0; i < Nvar; i++) {
			varint[i] = var0[i] + dtratio * (var1[i] - var0[i]);
		}
	}
}

void unscale_states(typeRNum *out, ctypeRNum *x, typeGRAMPC *grampc)
{
	char i;
	for (i = 0; i < grampc->Nx; i++) {
		out[i] = x[i] * grampc->xScale[i] + grampc->xOffset[i];
	}
}

void unscale_adjoints(typeRNum *out, ctypeRNum *adj, typeGRAMPC *grampc)
{
	char i;
	for (i = 0; i < grampc->Nx; i++) {
		out[i] = adj[i] / grampc->xScale[i];
	}
}

void unscale_controls(typeRNum *out, ctypeRNum *u, typeGRAMPC *grampc)
{
	char i;
	for (i = 0; i < grampc->Nu; i++) {
		out[i] = u[i] * grampc->uScale[i] + grampc->uOffset[i];
	}
}

void unscale_parameters(typeRNum *out, ctypeRNum *p, typeGRAMPC *grampc)
{
//	typeInt i;
/*	for (i = 0; i < grampc->Np; i++) {
		out[i] = p[i] * grampc->pScale[i] + grampc->pOffset[i];
	}*/
}

void unscale_time(typeRNum *out, ctypeRNum T, typeGRAMPC *grampc)
{
	out[0] = T * grampc->TScale + grampc->TOffset;
}


void scale_states(typeRNum *out, ctypeRNum *x, typeGRAMPC *grampc)
{
	char i;
	for (i = 0; i < grampc->Nx; i++) {
		out[i] = (x[i] - grampc->xOffset[i]) / grampc->xScale[i];
	}
}

void scale_adjoints(typeRNum *out, ctypeRNum *adj, typeGRAMPC *grampc)
{
	char i;
	for (i = 0; i < grampc->Nx; i++) {
		out[i] = adj[i] * grampc->xScale[i];
	}
}

void scale_controls(typeRNum *out, ctypeRNum *u, typeGRAMPC *grampc)
{
	char i;
	for (i = 0; i < grampc->Nu; i++) {
		out[i] = (u[i] - grampc->uOffset[i]) / grampc->uScale[i];
	}
}

void scale_parameters(typeRNum *out, ctypeRNum *p, typeGRAMPC *grampc)
{
//	typeInt i;
/*	for (i = 0; i < grampc->Np; i++) {
		out[i] = (p[i] - grampc->pOffset[i]) / grampc->pScale[i];
	}*/
}

void scale_time(typeRNum *out, typeRNum T, typeGRAMPC *grampc)
{
	out[0] = (T - grampc->TOffset) / grampc->TScale;
}

void scale_constraints(typeRNum *c, ctypeRNum *cScale, char Ncon)
{
	char i;
	for (i = 0; i < Ncon; i++) {
		c[i] = c[i] / cScale[i];
	}
}

void scale_cost(typeRNum *J, ctypeRNum JScale, ctypeInt Ncost)
{
	typeInt i;
	for (i = 0; i < Ncost; i++) {
		J[i] = J[i] / JScale;
	}
}
/*
void MatSetScalarp(typeRNum *C, ctypeRNum a, ctypeInt n1, ctypeInt n2)
{
	// matrix set C = a
	// C: (n1 x n2)
	typeInt i;
	for (i = 0; i < n1*n2; i++) {
		C[i] = 0.0;
	}
}
*/
void MatSetScalar(typeRNum *C, ctypeRNum a, ctypeInt n1, ctypeInt n2)
{
	/* matrix set C = a *
	 * C: (n1 x n2)     */
	typeInt i;
	for (i = 0; i < n1*n2; i++) {
		C[i] = a;
	}
}

void MatCopy(typeRNum *C, ctypeRNum *A, ctypeInt n1, ctypeInt n2)
{
	/* matrix copy C = A *
	 * A,C: (n1 x n2)    */
	//memcpy(C, A, n1*n2 * sizeof(*C));
	for (int i = 0; i<n1*n2; i++)
        C[i] = A[i];
}

void MatAdd(typeRNum *C, ctypeRNum *A, ctypeRNum *B, ctypeInt n1, ctypeInt n2)
{
	/* matrix summation C = A+B *
	 * A,B,C: (n1 x n2)         */
	typeInt i;
	for (i = 0; i < n1*n2; i++) {
		C[i] = A[i] + B[i];
	}
}

void MatMult(typeRNum *C, ctypeRNum *A, ctypeRNum *B, ctypeInt n1, ctypeInt n2, ctypeInt n3)
{
	/* matrix multiplication C = A*B			   *
	 * A: (n1 x n2)    B: (n2 x n3)   C: (n1 x n3) */
	typeInt i, j, k;
	typeRNum sigma;
	for (i = 0; i < n1; i++) {
		for (j = 0; j < n3; j++) {
			sigma = 0;
			for (k = 0; k < n2; k++) {
				sigma += A[i*n2 + k] * B[k*n3 + j];
			}
			C[i*n3 + j] = sigma;
		}
	}
}

void MatNorm(typeRNum *norm, ctypeRNum *A, ctypeInt n1, ctypeInt n2)
{
	/* A: (n1 x n2) */
	typeInt i, j;
	norm[0] = 0;
	for (i = 0; i < n1; i++) {
		for (j = 0; j < n2; j++) {
			norm[0] = norm[0] + A[i*n2 + j] * A[i*n2 + j];
		}
	}
	norm[0] = SQRT(norm[0]);
}

void MatDiffNorm(typeRNum *norm, ctypeRNum *A, ctypeRNum *B, ctypeInt n1, ctypeInt n2)
{
	/* A,B: (n1 x n2) */
	typeInt i, j;
	norm[0] = 0;
	for (i = 0; i < n1; i++) {
		for (j = 0; j < n2; j++) {
			norm[0] = norm[0] + (A[i*n2 + j] - B[i*n2 + j]) * (A[i*n2 + j] - B[i*n2 + j]);
		}
	}
	norm[0] = SQRT(norm[0]);
}
/*
typeInt grampc_estim_penmin(typeGRAMPC *grampc, ctypeInt rungrampc) {
	typeRNum Constraints[grampc->Nc], *xT_, *x_, *u_, *p_, h;
	typeRNum PenaltyMin_constr = grampc->PenaltyMin;
	typeRNum PenaltyMin_tol = grampc->PenaltyMin;
	typeRNum PenaltyMin = grampc->PenaltyMin;
	typeRNum NormTol = 0;
	typeRNum NormConstraints = 0;
	typeInt i, j, MaxGradIter, MaxMultIter, Status;

	// run grampc if flag is set to determine J and the trajectories
	if (rungrampc) {
		// limit the maximum number of iterations (important if e.g. the system is instable)
		MaxGradIter = grampc->MaxGradIter;
		MaxMultIter = grampc->MaxMultIter;
		if (MaxGradIter > 20) {
			grampc->MaxGradIter = 20;
		}
		if (MaxMultIter > 20) {
			grampc->MaxMultIter = 20;
		}
		grampc_run(grampc);
        }*/

	// assign parameters and final states
/*	if (grampc->ScaleProblem == INT_ON) {
		ASSIGN_P(p_, grampc->p, grampc);
		ASSIGN_X(xT_, grampc->x + grampc->Nx*(grampc->Nhor - 1), grampc);
        }
	else {
		p_ = grampc->p;

		xT_ = grampc->x + grampc->Nx*(grampc->Nhor - 1);
	//}
	//createNumMatrix(&Constraints, grampc->Nc);

	// evaluate and integrate the general constraints
 	for (i = 0; i < grampc->Nhor; i++) {
		if (grampc->ScaleProblem == INT_ON) {
			ASSIGN_X(x_, grampc->x + i * grampc->Nx, grampc);
			ASSIGN_U(u_, grampc->u + i * grampc->Nu, grampc);
		}
		else {
			x_ = grampc->x + i * grampc->Nx;
			u_ = grampc->u + i * grampc->Nu;
		}
		gfct(Constraints, grampc->t[i], x_, u_, p_, grampc->userparam);
		hfct(Constraints + grampc->Ng, grampc->t[i], x_, u_, p_, grampc->userparam);

		// scale constraints
		if (grampc->ScaleProblem == INT_ON) {
			scale_constraints(Constraints, grampc->cScale, grampc->Ng + grampc->Nh);
		}

		// integrate the squared constraints via trapezoidal rule
		if (i == 0) {
			h = (grampc->t[i + 1] - grampc->t[i]) / 2;
		}
		else if (i <= grampc->Nhor - 2) {
			h = (grampc->t[i + 1] - grampc->t[i - 1]) / 2;
		}
		else {
			h = (grampc->t[i] - grampc->t[i - 1]) / 2;
		}
		for (j = 0; j < grampc->Ng + grampc->Nh; j++) {
			NormConstraints += Constraints[j] * Constraints[j] * h;
		}
	}

	// evaluate the terminal constraints
	gTfct(Constraints + grampc->Ng + grampc->Nh, grampc->Thor, xT_, p_, grampc->userparam);
	hTfct(Constraints + grampc->Ng + grampc->Nh + grampc->NgT, grampc->Thor, xT_, p_, grampc->userparam);

	// scale constraints
 	if (grampc->ScaleProblem == INT_ON) {
		scale_constraints(Constraints + grampc->Ng + grampc->Nh, grampc->cScale + grampc->Ng + grampc->Nh, grampc->NgT + grampc->NhT);
	}

	// sum the squared terminal constarints
	for (i = grampc->Ng + grampc->Nh; i < grampc->Nc; i++) {
		NormConstraints += (Constraints[i] * Constraints[i]);
	}

	// set the constraints in relation to the costs
	if (NormConstraints > 0) {
		PenaltyMin_constr = 2 * ABS(grampc->J[0]) / NormConstraints;
	}
	//free(Constraints);

	// set the constraint tolerances in relation to the costs
 	for (i = 0; i < grampc->Nc; i++) {
		if (i < grampc->Ng + grampc->Nh) {
			NormTol = NormTol + grampc->Thor*grampc->ConstraintsAbsTol[i] * grampc->ConstraintsAbsTol[i];
		}
		else {
			NormTol = NormTol + grampc->ConstraintsAbsTol[i] * grampc->ConstraintsAbsTol[i];
		}
	}
	if (NormTol > 0) {
		PenaltyMin_tol = 2 * (typeRNum)1e-6 * ABS(grampc->J[0]) / NormTol;
	}

	// determine the resulting parameter and set the status
	if (PenaltyMin_constr > PenaltyMin_tol) {
		PenaltyMin = PenaltyMin_constr;
		Status = 1;
	}
	else {
		PenaltyMin = PenaltyMin_tol;
		Status = 2;
	}
	if (PenaltyMin > grampc->PenaltyMax / 500) {
		PenaltyMin = grampc->PenaltyMax / 500;
		Status = 0;
	}
	//grampc_setopt_real(grampc, "PenaltyMin", PenaltyMin);
	grampc->PenaltyMin = PenaltyMin;

	// reset grampc if grampc_run was evaluated
	if (rungrampc) {
		// reset iteration limits
		grampc->MaxGradIter = MaxGradIter;
		grampc->MaxMultIter = MaxMultIter;

		// reset sol structure
		MatSetScalar(grampc->xnext, 0, grampc->Nx, 1);
		MatSetScalar(grampc->unext, 0, grampc->Nu, 1);
//		MatSetScalar(grampc->pnext, 0, grampc->Np, 1);
		grampc->Tnext = 0;
		MatSetScalar(grampc->J, 0, 2, 1);
		grampc->cfct = 0;
		grampc->pen = 0;
		for (i = 0; i < grampc->MaxMultIter; i++) {
			grampc->iter[i] = 0;
		}
		grampc->status = 0;

		// reset rws structure
		MatSetScalar(grampc->tls, 0, grampc->Nhor, 1);

		MatSetScalar(grampc->x, 0, grampc->Nhor, grampc->Nx);
		MatSetScalar(grampc->adj, 0, grampc->Nhor, grampc->Nx);
		MatSetScalar(grampc->dcdx, 0, (grampc->Nhor + 1), grampc->Nx);

		MatSetScalar(grampc->uls, 0, grampc->Nhor, grampc->Nu);
		MatSetScalar(grampc->gradu, 0, grampc->Nhor, grampc->Nu);
		MatSetScalar(grampc->graduprev, 0, grampc->Nhor, grampc->Nu);
		MatSetScalar(grampc->dcdu, 0, grampc->Nhor, grampc->Nu);

		MatSetScalar(grampc->pls, 0, grampc->Np, 1);
		MatSetScalar(grampc->gradp, 0, grampc->Np, 1);
		MatSetScalar(grampc->gradpprev, 0, grampc->Np, 1);
		MatSetScalar(grampc->dcdp, 0, (grampc->Nhor + 1), grampc->Np);

		grampc->gradT = 0;
		grampc->gradTprev = 0;
		grampc->dcdt = 0;

		MatSetScalar(grampc->mult, 0, grampc->Nhor, grampc->Nc);
		MatSetScalar(grampc->rws_pen, 0, grampc->Nhor, grampc->Nc);
		MatSetScalar(grampc->rws_cfct, 0, grampc->Nhor, grampc->Nc);
		MatSetScalar(grampc->cfctprev, 0, grampc->Nhor, grampc->Nc);

		MatSetScalar(grampc->rwsScale, 0, 2 * (grampc->Nx + grampc->Nu + grampc->Np), 1);

		// Initialization
		init_rws_time(&grampc);
		init_rws_controls(&grampc);
//		init_rws_parameters(grampc);
		init_rws_multipliers(&grampc);
		init_rws_linesearch(&grampc);

//		resize_rwsRodas(grampc);
		resize_rwsGeneral(&grampc);
	}
	return Status;
}*/


#ifdef __cplusplus
}
#endif // __cplusplus
