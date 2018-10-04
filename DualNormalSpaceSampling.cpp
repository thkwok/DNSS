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

#include <algorithm>

#define M_PI 3.1415926535897932384626433832795

#define UseAverage 0

int Nsize[2] = { 12,6 };
int Rsize[2] = { 6,6 };

void VectorProduct(double n1[], double n2[], double n3[])
{
	n3[0] = n1[1] * n2[2] - n1[2] * n2[1];
	n3[1] = n1[2] * n2[0] - n1[0] * n2[2];
	n3[2] = n1[0] * n2[1] - n1[1] * n2[0];
}

double VectorProject(double n1[], double n2[])
{
	double r = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
	return r;
}

double CalLength(double n[], int d=3)
{
	double tt = 0;
	for (int i = 0; i<d; i++) tt += n[i] * n[i];
	tt = sqrt(tt);
	return tt;
}

bool Normalize(double n[], int d=3)
{
	double tt = CalLength(n, d);
	if (tt<1e-5) {
		for (int i = 0; i<d; i++) n[i] = 0;
		return false;
	}
	else if (fabs(tt - 1)<1e-5) {
		return true;
	}
	else {
		for (int i = 0; i<d; i++) n[i] /= tt;
	}
	return true;
}


class GLKHeapNode
{
public:
	GLKHeapNode() { attachedObj = NULL; };
	~GLKHeapNode() {};
	double GetValue() { return v; };
	void SetValue(double value) { v = value; };

	void* attachedObj;
	int i;
	double v;
};

bool cmpLess(GLKHeapNode *a, GLKHeapNode *b) {
	return a->GetValue() > b->GetValue();
}



DualNormalSpaceSampling::DualNormalSpaceSampling(list<MeshNode>& nList)
{
	nodeList = &nList;
	Bin = 0; BinR = 0;

	int i = 0;
	for (std::list<MeshNode>::iterator it = nList.begin(); it != nList.end(); ++it) {
		it->index = i++;
	}
}

DualNormalSpaceSampling::~DualNormalSpaceSampling(void)
{
	ClearBin();
}


void DualNormalSpaceSampling::ClearBin()
{
	if (!Bin) return;
	for (int i = 0; i < Nsize[0]; i++) delete[] Bin[i]; delete[] Bin; Bin = 0;
	for (int i = 0; i < Rsize[0]; i++) delete[] BinR[i]; delete[] BinR; BinR = 0;
}

void DualNormalSpaceSampling::InitialBin()
{
	ClearBin();

	Bin = new list<GLKHeapNode*>*[Nsize[0]];
	for (int i = 0; i < Nsize[0]; i++) {
		Bin[i] = new list<GLKHeapNode*>[Nsize[1]];
	}

	BinR = new list<GLKHeapNode*>*[Rsize[0]];
	for (int i = 0; i < Rsize[0]; i++) {
		BinR[i] = new list<GLKHeapNode*>[Rsize[1]];
	}

	heapnodeArr.resize(nodeList->size());
}


void DualNormalSpaceSampling::FindBinNumber(double normal[3], int binNum[2])
{
	//compute theta
	double theta = 0;
	double x = normal[0], y = normal[1];
	double l = sqrt(x*x + y*y);
	if (l != 0) {
		x /= l; y /= l;
		theta = acos(fabs(y));
		if (x >= 0 && y >= 0) {}
		else if (x >= 0 && y < 0) theta = M_PI - theta;
		else if (x < 0 && y < 0) theta += M_PI;
		else if (x < 0 && y >= 0) theta = 2 * M_PI - theta;
	}

	double z = normal[2];

	//convert to u and v
	binNum[0] = (int)(theta / (2* M_PI) * Nsize[0]);
	binNum[1] = (int)((z + 1) / 2 * Nsize[1]);
	
	if (binNum[0] == Nsize[0]) binNum[0]--;
	if (binNum[1] == Nsize[1]) binNum[1]--;
}

void DualNormalSpaceSampling::FindBinNumberRot(double center[3], double p[3], double normal[3], int binNum[2])
{
	double theta = 0;
	double op[3]; for (int i = 0; i < 3; i++) op[i] = p[i] - center[i];
	double r[3]; VectorProduct(op, normal, r);
	Normalize(r);
	double x = r[0], y = r[1], z = r[2];
	double l = sqrt(x*x + y*y); x /= l; y /= l;
	theta = acos(fabs(y));
	if (x >= 0 && y >= 0) {}
	else if (x >= 0 && y < 0) theta = M_PI - theta;
	else if (x < 0 && y < 0) {}
	else if (x < 0 && y >= 0) theta = M_PI - theta;

	binNum[0] = (int)(theta / M_PI * Rsize[0]);
	binNum[1] = (int)((z + 1) / 2 * Rsize[1]);

	if (binNum[0] == Nsize[0]) binNum[0]--;
	if (binNum[1] == Nsize[1]) binNum[1]--;
}

double DualNormalSpaceSampling::ComputeMu(double center[3], double p[3], double n[3], double r[3], double normalization)
{
	double theta = M_PI / 4.0;

	// get beta
	double op[3]; for (int i = 0; i < 3; i++) op[i] = p[i] - center[i];
	double lop = CalLength(op);

	double dot = VectorProject(op, n) / lop;
	double beta = acos(dot);
	if (beta < 0) beta += M_PI;

	// calculate rotational return
	double pp = sqrt(2 * (1 - cos(theta)));

	double a[2] = { beta, M_PI - beta };
	double m[2] = { 0 };

	for (int i = 0; i < 2; i++) {
		double pq = pp * cos(a[i] - 0.5*theta);
		double hp = pq * sin(a[i]);
		double hq = pq * cos(a[i]);

		m[i] = (tan(theta) - hq / (1 - hp)) / tan(theta) * lop / normalization;
	}

	double mu = (m[0] > m[1]) ? m[0] : m[1];
	return mu;
}

void DualNormalSpaceSampling::InsertNodeToList(list<GLKHeapNode*> &list, GLKHeapNode *heapnode)
{
	double mu = heapnode->GetValue();

	if (list.empty()) {
		list.push_back(heapnode); return;
	}

	//quick check the end node
	{
		GLKHeapNode *node = (GLKHeapNode*)list.back();
		double v = node->GetValue();
		if (mu <= v) {
			list.push_back(heapnode);
			return;
		}
	}

	for (std::list<GLKHeapNode*>::iterator it = list.begin(); it != list.end(); ++it) {
		GLKHeapNode *node = (GLKHeapNode*)*it;
		double v = node->GetValue();
		if (mu >= v) {
			if (it == list.begin()) list.push_front(heapnode);
			else list.insert(++it, heapnode);
			return;
		}
	}

	printf("ERROR: Cannot add heapnode?\n");
}

int DualNormalSpaceSampling::AssignPointsToBin()
{
	InitialBin();

	int NumPts = nodeList->size();
	double center[3] = { 0 };
	for (auto const& node : *nodeList){
		for (int i = 0; i < 3; i++) center[i] += node.v[i];
	}
	for (int i = 0; i < 3; i++) center[i] /= double(NumPts);
	//printf("center: %f %f %f\n", center[0], center[1], center[2]);

	double Lmax = 0, Lavg = 0; //for normalization
	for (auto const& node : *nodeList) {
		double p[3]; for (int i = 0; i < 3; i++) p[i] = node.v[i] - center[i];
		double l = CalLength(p);
		if (l > Lmax) Lmax = l;
		Lavg += l;
	}
	Lavg /= NumPts;
	if (UseAverage) Lmax = Lavg;

	int count = 0; int id = -1;
	for (auto const& node : *nodeList) {
		++id;
		double p[3], n[3];
		for (int i = 0; i < 3; i++) {
			p[i] = node.v[i]; n[i] = node.n[i];
		}
		double op[3]; for (int i = 0; i < 3; i++) op[i] = p[i] - center[i];
		double r[3]; VectorProduct(op, n, r);
		Normalize(r);

		double lop = CalLength(op);

		double mu = ComputeMu(center, p, n, r, Lmax);
		if (mu < 0) continue;

		heapnodeArr[id] = new GLKHeapNode();
		GLKHeapNode *heapnode = heapnodeArr[id];
		heapnode->attachedObj = (void*)&node;
		heapnode->SetValue(mu); 
		heapnode->i = 0; // means not selected

		int b[2]; FindBinNumber(n, b);
		InsertNodeToList(Bin[b[0]][b[1]], heapnode);

		int bR[2]; FindBinNumberRot(center, p, n, bR);
		InsertNodeToList(BinR[bR[0]][bR[1]], heapnode);

		count++;
	}
	return count;
}

bool DualNormalSpaceSampling::SamplingOnRotSpace(list<MeshNode> *outputList, int NumSample)
{
	for (int i = 0; i < Rsize[0]; i++) for (int j = 0; j < Rsize[1]; j++) {
		list<GLKHeapNode*> &list = BinR[i][j];

		if (list.empty()) continue;
		GLKHeapNode* heapnode = list.front(); list.pop_front();
		while (heapnode->i) {
			if (list.empty()) { heapnode = 0; break; }
			heapnode = list.front(); list.pop_front();
		} if (!heapnode) continue;


		MeshNode &node = *(MeshNode*)heapnode->attachedObj;
		double mu = heapnode->GetValue();
		heapnode->i = 1;

		outputList->push_back(node);
		//printf("added one node from bin %d %d %d\n", b[0], b[1], b[2]);

		if (outputList->size() >= NumSample) { i = Rsize[0]; j = Rsize[1]; break; }
	}

	return (outputList->size() == NumSample);
}

bool DualNormalSpaceSampling::SamplingOnNormalSpace(list<MeshNode> *outputList, int NumSample)
{
	for (int i = 0; i < Nsize[0]; i++) for (int j = 0; j < Nsize[1]; j++) {
		list<GLKHeapNode*> &list = Bin[i][j];

		if (list.empty()) continue;
		GLKHeapNode* heapnode = list.front(); list.pop_front();
		while (heapnode->i) {
			if (list.empty()) { heapnode = 0; break; }
			heapnode = list.front(); list.pop_front();
		} if (!heapnode) continue;

		MeshNode &node = *(MeshNode*)heapnode->attachedObj;
		double mu = heapnode->GetValue();
		heapnode->i = 1;

		outputList->push_back(node);
		//printf("added one node from bin %d %d %d\n", b[0], b[1], b[2]);

		if (outputList->size() >= NumSample) { i = Nsize[0]; j = Nsize[1]; break; }
	}

	return (outputList->size() == NumSample);
}

bool DualNormalSpaceSampling::SamplingOnScrewSpace(list<MeshNode> *outputList, int NumSample)
{
	//first pass on Rot Space
	SamplingOnRotSpace(outputList, NumSample);
	if (outputList->size() == NumSample) return true;

	//if first pass isn't full, then use histogram to guide the point picking

	//get center
	int NumPts = nodeList->size();
	double center[3] = { 0 };
	for (auto const &node : *nodeList){
		for (int i = 0; i < 3; i++) center[i] += node.v[i];
	}
	for (int i = 0; i < 3; i++) center[i] /= double(NumPts);
	//printf("center: %f %f %f\n", center[0], center[1], center[2]);


	//build histogram
	double **H = new double*[Nsize[0]];
	for (int i = 0; i < Nsize[0]; i++) { H[i] = new double[Nsize[1]]; for (int j = 0; j < Nsize[1]; j++) H[i][j] = 0;}
	double **HR = new double*[Rsize[0]];
	for (int i = 0; i < Rsize[0]; i++) { HR[i] = new double[Rsize[1]]; for (int j = 0; j < Rsize[1]; j++) HR[i][j] = 0;}

	BuildHistogram(outputList, H, HR);
	
	//make two heap to guide the process
	int NumT = Nsize[0] * Nsize[1], NumR = Rsize[0]*Rsize[1];
	vector<GLKHeapNode*> heap, heapR;
	heap.reserve(NumT); heapR.reserve(NumR);
	GLKHeapNode **hnode = new GLKHeapNode*[NumT], **hnodeR = new GLKHeapNode*[NumR];
	for (int i = 0; i < NumT; i++) {
		hnode[i] = 0;
		int x = i / Nsize[1], y = i % Nsize[1];
		list<GLKHeapNode*> &list = Bin[x][y];
		if (list.empty()) continue; //don't put in heap if there are no points can be added in that bin
		hnode[i] = new GLKHeapNode();
		hnode[i]->SetValue(H[x][y]); hnode[i]->i = i;
		hnode[i]->attachedObj = (void*)&list;
		heap.push_back(hnode[i]);
	}
	for (int i = 0; i < NumR; i++) {
		hnodeR[i] = 0;
		int x = i / Rsize[1], y = i % Rsize[1];
		list<GLKHeapNode*> &list = BinR[x][y];
		if (list.empty()) continue; //don't put in heap if there are no points can be added in that bin
		hnodeR[i] = new GLKHeapNode();
		hnodeR[i]->SetValue(HR[x][y]); hnodeR[i]->i = i;
		hnodeR[i]->attachedObj = (void*)&list;
		heapR.push_back(hnodeR[i]);
	}

	make_heap(heap.begin(), heap.end(), cmpLess); 
	make_heap(heapR.begin(), heapR.end(), cmpLess);

	//pick points
	while (outputList->size() < NumSample) {
		//check the minimum is from which heap
		if (heapR.front()->GetValue() > heap.front()->GetValue()) {  
			//for translation
			GLKHeapNode *hn = heap.front();
			int id = hn->i;
			int x = id / Nsize[1], y = id % Nsize[1];
			//pick the lowest one to add point
			list<GLKHeapNode*> &hlist = *(list<GLKHeapNode*>*)hn->attachedObj;

			if (hlist.empty()) { pop_heap(heap.begin(), heap.end(), cmpLess); heap.pop_back(); continue; }
			GLKHeapNode* heapnode = hlist.front(); hlist.pop_front();
			while (heapnode->i) {
				if (hlist.empty()) { heapnode = 0; break; }
				heapnode = hlist.front(); hlist.pop_front();
			} if (!heapnode) { pop_heap(heap.begin(), heap.end(), cmpLess); heap.pop_back(); continue; }

			MeshNode &node = *(MeshNode*)heapnode->attachedObj;
			double mu = heapnode->GetValue();
			heapnode->i = 1;

			outputList->push_back(node);

			//update histograms and heaps
			H[x][y]+=1;
			if (hlist.empty()) { pop_heap(heap.begin(), heap.end(), cmpLess); heap.pop_back(); }
			else { 
				pop_heap(heap.begin(), heap.end(), cmpLess);
				hn->SetValue(H[x][y]);
				push_heap(heap.begin(), heap.end(), cmpLess);
			}
			//rotation by the node
			double n[3], p[3];
			for (int i = 0; i < 3; i++) {
				n[i] = node.n[i]; p[i] = node.v[i];
			}
			int bR[2]; FindBinNumberRot(center, p, n, bR);
			HR[bR[0]][bR[1]] += mu;
			GLKHeapNode *updateNode = hnodeR[bR[0] * Rsize[1] + bR[1]];
			pop_heap(heapR.begin(), heapR.end(), cmpLess);
			updateNode->SetValue(HR[bR[0]][bR[1]]);
			push_heap(heapR.begin(), heapR.end(), cmpLess);
		}
		else { //for rotation
			GLKHeapNode *hn = heapR.front();
			int id = hn->i;
			int x = id / Rsize[1], y = id % Rsize[1];
			//pick the lowest one to add point
			list<GLKHeapNode*> &hlist = *(list<GLKHeapNode*>*)hn->attachedObj;

			if (hlist.empty()) { pop_heap(heapR.begin(), heapR.end(), cmpLess); heapR.pop_back(); continue; }
			GLKHeapNode* heapnode = hlist.front(); hlist.pop_front();
			while (heapnode->i) {
				if (hlist.empty()) { heapnode = 0; break; }
				heapnode = hlist.front(); hlist.pop_front();
			} if (!heapnode) { pop_heap(heapR.begin(), heapR.end(), cmpLess); heapR.pop_back(); continue; }

			MeshNode &node = *(MeshNode*)heapnode->attachedObj;
			double mu = heapnode->GetValue();
			heapnode->i = 1;

			outputList->push_back(node);

			//update histograms and heaps
			HR[x][y]+= mu;
			if (hlist.empty()) { pop_heap(heapR.begin(), heapR.end(), cmpLess); heapR.pop_back(); }
			else { 
				pop_heap(heapR.begin(), heapR.end(), cmpLess);
				hn->SetValue(HR[x][y]); 
				push_heap(heapR.begin(), heapR.end(), cmpLess);
			}
			//translation by the node
			double n[3]; for (int i = 0; i < 3; i++) n[i] = node.n[i];
			int b[2]; FindBinNumber(n, b);
			H[b[0]][b[1]]+=1;
			GLKHeapNode *updateNode = hnode[b[0] * Nsize[1] + b[1]];
			pop_heap(heap.begin(), heap.end(), cmpLess);
			updateNode->SetValue(H[b[0]][b[1]]);
			push_heap(heap.begin(), heap.end(), cmpLess);
		}
	}


	delete[]hnode; 
	delete[]hnodeR;
	heap.clear();
	heapR.clear();

	/////////////////////////////////////
	for (int i = 0; i < Nsize[0]; i++) delete[] H[i]; delete[]H;
	for (int i = 0; i < Rsize[0]; i++) delete[] HR[i]; delete[]HR;

	return (outputList->size() == NumSample);
}

bool DualNormalSpaceSampling::Run(list<int> &SampleIndex, int NumSample)
{
	SampleIndex.clear();
	if (nodeList->empty()) return false;
	if (NumSample > nodeList->size()){
		for (int i = 0; i < nodeList->size(); i++) SampleIndex.push_back(i);
		return true;
	}

	int count = AssignPointsToBin();

	if (count < NumSample) {
		printf("There is not enough points for samples needed. (%d\%d) \n", count, NumSample);
		NumSample = count;
	}

	list<MeshNode> outputList;

	SamplingOnScrewSpace(&outputList, NumSample);

	for (auto const &node : outputList) {
		SampleIndex.push_back(node.index);
	}

	return (outputList.size() == NumSample);
}


void DualNormalSpaceSampling::BuildHistogram(list<MeshNode> *list, double **H, double **HR)
{
	//get center
	int NumPts = nodeList->size();
	double center[3] = { 0 };
	for (auto const &node : *nodeList){
		for (int i = 0; i < 3; i++) center[i] += node.v[i];
	}
	for (int i = 0; i < 3; i++) center[i] /= double(NumPts);
	//printf("center: %f %f %f\n", center[0], center[1], center[2]);

	double Lmax = 0; //for normalization
	double Lavg = 0;
	for (auto const &node : *nodeList) {
		double p[3];
		for (int i = 0; i < 3; i++) p[i] = node.v[i]-center[i];
		double l = CalLength(p);
		if (l > Lmax) Lmax = l;
		Lavg += l;
	} 
	Lavg /= NumPts;
	if (UseAverage) Lmax = Lavg;
	//printf("Lmax = %f\n", Lmax);

	for (auto const &node : *list){
		double n[3], p[3];
		for (int i = 0; i < 3; i++) {
			n[i] = node.n[i]; p[i] = node.v[i];
		}

		double op[3]; for (int i = 0; i < 3; i++) op[i] = p[i] - center[i];
		double r[3]; VectorProduct(op, n, r);
		Normalize(r);

		double mu = ComputeMu(center, p, n, r, Lmax);
		//if (mu < 0) continue;

		int b[2]; FindBinNumber(n, b);
		int bR[2]; FindBinNumberRot(center, p, n, bR);

		H[b[0]][b[1]]+=1.0;
		HR[bR[0]][bR[1]] += mu;
	}

}