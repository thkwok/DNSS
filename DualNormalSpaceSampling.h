/*
 *  Copyright (C) 2018, Tsz-Ho Kwok in Concordia University, Montreal
 *  All rights reserved.
 *  
 *   
 *  Redistribution and use in source and binary forms, with or without modification, 
 *  are permitted provided that the following conditions are met:
 *  
 *  1. Redistributions of source code must retain the above copyright notice, 
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice, 
 *     this list of conditions and the following disclaimer in the documentation 
 *	   and/or other materials provided with the distribution.
 *  3. Acknowledge the provided code by citing the paper in any publication using this code.
 *     T.-H. Kwok, "DNSS: Dual-Normal Space Sampling for 3D ICP Registration", IEEE TASE, 2018.
 *
 *   THIS CODE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
 *   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
 *   IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
 *   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
 *   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
 *   OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
 *   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY 
 *   OF SUCH DAMAGE.
 */


#pragma once

#include <list>
#include <vector>       // std::vector
using namespace std;

class GLKHeapNode;

///////////////////////////////////////////////////
// define node structure -- with vertex position and normal, as well as its index
struct MeshNode {
	double v[3]; double n[3]; int index;
};
///////////////////////////////////////////////////

class DualNormalSpaceSampling
{
public:
	DualNormalSpaceSampling(list<MeshNode>& nList);
	~DualNormalSpaceSampling(void);

	bool Run(list<int> &SampleIndex, int NumSample);

private:
	int AssignPointsToBin();


	void ClearBin();
	void InitialBin(); //initial bin for normal space (theta is 12, phi is 6)
	void FindBinNumber(double normal[3], int binNum[2]);
	void FindBinNumberRot(double center[3], double p[3], double normal[3], int binNum[2]);

	double ComputeMu(double center[3], double p[3], double n[3], double r[3], double normalization = 1);

	void InsertNodeToList(list<GLKHeapNode*> &list, GLKHeapNode *heapnode);

	void BuildHistogram(list<MeshNode> *list, double **H, double **HR);


	bool SamplingOnRotSpace(list<MeshNode> *outputList, int NumSample);
	bool SamplingOnNormalSpace(list<MeshNode> *outputList, int NumSample);
	bool SamplingOnScrewSpace(list<MeshNode> *outputList, int NumSample);

private:
	list<MeshNode> *nodeList;

	list<GLKHeapNode*>** Bin;
	list<GLKHeapNode*>** BinR;

	vector <GLKHeapNode*> heapnodeArr;
};
