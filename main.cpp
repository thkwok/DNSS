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

#include "DualNormalSpaceSampling.h"

int main()
{
	double x, y, z; //point positions
	double nx, ny, nz; //point normals
	list<MeshNode> nodeList;
	//store points by MeshNode
	for (...) {
		MeshNode node;
		node.v[0] = x; node.v[1] = y; node.v[2] = z;
		node.n[0] = nx; node.n[1] = ny; node.n[2] = nz;
		nodeList.push_back(m_node);
	}

	DualNormalSpaceSampling dnss(nodeList);

	list<int> SampleIndex; //the output -- stored the index of the nodes (start from 0)

	int NumPts = 0; //to get the target number of sample points
	printf("Number of sampling points: "); scanf("%d", &NumPts);

	if (dnss.Run(SampleIndex, NumPts)) {
		printf("Selected %d points in Dual Normal Space Sampling\n", SampleIndex.size());
		printf("the indices of selected points are: ");
		for (auto const &index : SampleIndex) {
			printf("index ");
		}printf("\n");
	}
	else printf("ERROR: Cannot select points by Dual Normal Space Sampling1\n");


	return 0;
}