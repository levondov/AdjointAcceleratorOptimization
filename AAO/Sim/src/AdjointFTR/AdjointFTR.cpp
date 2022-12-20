// AdjointFTR.cpp : Defines the entry point for the application.
//

#include "AdjointFTR.h"
#include <math.h>

using namespace std;

namespace FTR
{
	AdjointFTR::AdjointFTR()
	{

	}

	AdjointFTR::~AdjointFTR()
	{

	}

	double* AdjointFTR::GetONmats(double kPerv, double kSol, double kQuad, double kQuadRot, double pipeRadius, double Y0, double Y1, double Y2, double Y10)
	{
		double cq = cos(2.0 * Y10 - 2.0 * kQuadRot);
		double sq = sin(2.0 * Y10 - 2.0 * kQuadRot);

		// calculate predefined variables
		double dotTerm = Y0 * Y0 - Y1 * Y1 - Y2 * Y2;
		if (dotTerm <= 0.0) { dotTerm = 1.0e-16; }
		double qDelta = sqrt(dotTerm);
		double ab4 = 1.0 / qDelta;
		double ca_ab4 = -Y1 / ((Y0 + qDelta) * qDelta);
		double sa_ab4 = -Y2 / ((Y0 + qDelta) * qDelta);

		// pipe constant
		double pipeConstant;
		double omatPipeConstantTerm[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		if (pipeRadius != 0)
		{
			pipeConstant = (8.0 * kPerv) / (pipeRadius * pipeRadius * pipeRadius * pipeRadius);
			omatPipeConstantTerm[1] = -Y1 * pipeConstant;
			omatPipeConstantTerm[2] = -Y2 * pipeConstant;
			omatPipeConstantTerm[3] = -Y1 * pipeConstant;
			omatPipeConstantTerm[6] = -Y2 * pipeConstant;
			omatPipeConstantTerm[10] = -Y2 * pipeConstant;
			omatPipeConstantTerm[11] = Y1 * pipeConstant;
		}
		else
		{
			pipeConstant = 0;
		}

		double* onmats = new double[12];
		onmats[0] = -1 * (kSol * kSol) / 2.0 + ab4 * kPerv + omatPipeConstantTerm[0];
		onmats[1] = 2.0 * kQuad * cq + ca_ab4 * kPerv + omatPipeConstantTerm[1];
		onmats[2] = -2.0 * kQuad * sq + sa_ab4 * kPerv + omatPipeConstantTerm[2];
		onmats[3] = 2.0 * kQuad * cq + ca_ab4 * kPerv + omatPipeConstantTerm[3];
		onmats[4] = -1 * (kSol * kSol) / 2.0 + ab4 * kPerv + omatPipeConstantTerm[4];
		onmats[5] = 0.0 + omatPipeConstantTerm[5];
		onmats[6] = -2.0 * kQuad * sq + sa_ab4 * kPerv + omatPipeConstantTerm[6];
		onmats[7] = 0.0 + omatPipeConstantTerm[7];
		onmats[8] = -1 * (kSol * kSol) / 2.0 + ab4 * kPerv + omatPipeConstantTerm[8];

		onmats[9] = 0.0 - omatPipeConstantTerm[9];
		onmats[10] = 2.0 * kQuad * sq - sa_ab4 * kPerv - omatPipeConstantTerm[10];
		onmats[11] = 2.0 * kQuad * cq + ca_ab4 * kPerv - omatPipeConstantTerm[11];

		return onmats;
	}

	double* AdjointFTR::GetSCVM(double kPerv, double Y0, double Y1, double Y2, double Y3, double Y4, double Y5)
	{
		double* svcm = new double[27];
		int svcmCC = 0;

		// calculate predefined variables
		double dotTerm = Y0 * Y0 - Y1 * Y1 - Y2 * Y2;
		if (dotTerm <= 0.0) { dotTerm = 1.0e-16; }
		double qDelta = sqrt(dotTerm);
		double qDeltaPlus = qDelta + Y0;

		// V vectors
		double V1[3] = {}; double V2[3] = {}; double V3[3] = {}; double V4[3] = {};
		double V1_c = -1 * (kPerv / (qDelta * qDelta));
		double V2_c = (kPerv / (qDelta * qDeltaPlus * qDeltaPlus));
		double V3_c = -1 * (kPerv / (qDelta * qDeltaPlus));
		V1[0] = V1_c * (Y3 - Y1 * Y4 / qDeltaPlus - Y2 * Y5 / qDeltaPlus);
		V1[1] = V1_c * (-Y1 * Y3 / qDeltaPlus + Y4);
		V1[2] = V1_c * (-Y2 * Y3 / qDeltaPlus + Y5);
		V2[0] = V2_c * (Y1 * Y4 + Y2 * Y5);
		V2[1] = V2_c * (Y1 * Y3);
		V2[2] = V2_c * (Y2 * Y3);
		V3[0] = V3_c * Y4;
		V3[1] = V3_c * Y3;
		V3[2] = 0.0;
		V4[0] = V3_c * Y5; //V4_c = V3_c
		V4[1] = 0.0;
		V4[2] = V3_c * Y3; //V4_c = V3_c

		// U vectors
		double U1[3] = {}; double U2[3] = {}; double U3[3] = {}; double U4[3] = {};
		double U1_c = 1 / qDelta;
		U1[0] = U1_c * Y0;
		U1[1] = U1_c * -Y1;
		U1[2] = U1_c * Y2;
		U2[0] = U1[0] + 1.0;
		U2[1] = U1[1];
		U2[2] = U1[2];
		U3[0] = 0.0;
		U3[1] = 1.0;
		U3[2] = 0.0;
		U4[0] = 0.0;
		U4[1] = 0.0;
		U4[2] = 1.0;

		// W vectors
		double W1[3] = {}; double W2[3] = {}; double W3[3] = {}; double W4[3] = {};
		double W1_c = -1 * (kPerv / qDelta);
		W1[0] = W1_c * 1.0;
		W1[1] = W1_c * (Y1 / qDeltaPlus);
		W1[2] = W1_c * (Y2 / qDeltaPlus);
		W2[0] = V2_c * (Y1 * Y1 + Y2 * Y2); // W2_c = V2_c
		W2[1] = V2_c * (Y1 * Y0);			// W2_c = V2_c
		W2[2] = V2_c * (Y2 * Y0);			// W2_c = V2_c
		W3[0] = V3_c * Y1;					// W3_c = V3_c
		W3[1] = V3_c * Y0;					// W3_c = V3_c 
		W3[2] = 0.0;
		W4[0] = V3_c * Y2;					// W4_c = V3_c 
		W4[1] = 0.0;						// W4_c = V3_c 
		W4[2] = V3_c * Y0;					// W4_c = V3_c 

		// X vectors
		double X1[3] = {}; double X2[3] = {}; double X3[3] = {}; double X4[3] = {};
		X1[0] = 0.0;
		X1[1] = -1 * kPerv * Y2 / (qDeltaPlus * qDelta * qDelta);
		X1[2] = kPerv * Y1 / (qDeltaPlus * qDelta * qDelta);
		X2[0] = 0.0;
		X2[1] = -1 * kPerv * Y2 / (qDeltaPlus * qDeltaPlus * qDelta);
		X2[2] = kPerv * Y1 / (qDeltaPlus * qDeltaPlus * qDelta);
		X3[0] = 0.0;
		X3[1] = 0.0;
		X3[2] = -1 * kPerv / (qDeltaPlus * qDelta);
		X4[0] = 0.0;
		X4[1] = kPerv / (qDeltaPlus * qDelta);
		X4[2] = 0.0;

		// outer matrix multiply
		// Mq = W1*U1 + W2*U2 + W3*U3 + W4*U4
		// Mp = V1*U1 + V2*U2 + V3*U3 + V4*U4
		// Mn = X1*U1 + X2*U2 + X3*U3 + X4*U4
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				svcm[svcmCC] = W1[i] * U1[j] + W2[i] * U2[j] + W3[i] * U3[j] + W4[i] * U4[j];
				svcmCC++;
			}
		}
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				svcm[svcmCC] = V1[i] * U1[j] + V2[i] * U2[j] + V3[i] * U3[j] + V4[i] * U4[j];
				svcmCC++;
			}
		}
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				svcm[svcmCC] = X1[i] * U1[j] + X2[i] * U2[j] + X3[i] * U3[j] + X4[i] * U4[j];
				svcmCC++;
			}
		}

		return svcm;
	}
}

extern "C"
{
	__declspec(dllexport) FTR::AdjointFTR* AdjointFTR_new()
	{
		return new FTR::AdjointFTR();
	}
	__declspec(dllexport) void AdjointFTR_destroy(FTR::AdjointFTR* adjointFTR)
	{
		delete adjointFTR;
	}
	__declspec(dllexport) double* AdjointFTR_getSCVM(FTR::AdjointFTR* adjointFTR, double kPerv, double Y0, double Y1, double Y2, double Y3, double Y4, double Y5)
	{
		return adjointFTR->GetSCVM(kPerv, Y0, Y1, Y2, Y3, Y4, Y5);
	}
	__declspec(dllexport) double* AdjointFTR_getONmats(FTR::AdjointFTR* adjointFTR, double kPerv, double kSol, double kQuad, double kQuadRot, double pipeRadius, double Y0, double Y1, double Y2, double Y10)
	{
		return adjointFTR->GetONmats(kPerv, kSol, kQuad, kQuadRot, pipeRadius, Y0, Y1, Y2, Y10);
	}
}
