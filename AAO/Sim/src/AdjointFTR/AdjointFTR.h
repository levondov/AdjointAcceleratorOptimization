// AdjointFTR.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <iostream>

#ifndef AdjointFTR_H
#define AdjointFTR_H

namespace FTR
{
	class AdjointFTR
	{

	public:

		AdjointFTR();
		~AdjointFTR();

		double* GetSCVM(double kPerv, double Y0, double Y1, double Y2, double Y3, double Y4, double Y5);
		double* GetONmats(double kPerv, double kSol, double kQuad, double kQuadRot, double pipeRadius, double Y0, double Y1, double Y2, double Y10);

	private:

	};
}

#endif