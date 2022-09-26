#include "cplex.h"
#include <iomanip>
#include <iostream>
#include <string> 
#include <fstream>
#include <sstream>  //istringstream必须包含这个头文件
#include <algorithm>
#include <math.h>
#include <list>
#include <vector> 
#include <cfloat>
#include <ctime>
#include<cstdlib>
#include <limits.h>

using namespace std;

void reset_all();
void read_data(string fdevice, string faccess, string fSNR);
void init(string fdevice, string faccess, string fSNR);
void calculate_price(char p, int i = 0, int k = 0);
void init_UL();
bool is_full_cpubw(int i, int j);
bool is_full_bw(int i, int j);
bool is_full(int j);
void cal_energy();
void EB_CMA(int j);
void print_result(Result* cplex_r, Result* cplexcloud_r, Result* cplexedge_r, Result* cma_r, Result* ff_r, Result* e_r);
void print_csv(Result* cplex_r, Result* cplexcloud_r, Result* cplexedge_r, Result* cma_r, Result* ff_r, Result* e_r);
void print_data();
void first_fit_alg(Result* result);
void equal_alg(Result* result);
void redefinition();

// 在实验每组数据之后，重置全局变量
void reset_all()
{
	for (int i = 0; i < ACCESS; i++) {
		access[i].constrain_BW = 0;
		access[i].sigma = 0;
		access[i].v = 0;
		access[i].p = 0;
		access[i].z_BW = 0;
		for (int k = 0; k < RESOURCE; k++) {
			access[i].constrain_r[k] = 0;
			access[i].z_r[k] = 0;
			U_a[i][k] = 0.0f;
			L_a[i][k] = (double)INF;
			alpha[i][k] = 0.0f;
		}
		for (int j = 0; j < DEVICES; j++) {
			x[i][j] = 0;
			y[i][j] = 0;
			f_edge[i][j] = 0; //i为j分配的算力
			time_tran[i][j] = 0;
			energy_tran_access[i][j] = 0;
			energy_comp_edge[i][j] = 0;
			energy_tran_cloud1[i][j] = 0;
			Energy_edge[i][j] = 0;
			Energy_cloud[i][j] = 0;
			eta_edge[i][j] = 0;
			eta_cloud[i][j] = 0;

			energy_edgeserver1[i][j] = 0;  
			energy_tran1[i][j] = 0;
		}
		energy_tran_cloud2[i] = 0;
		U_a_BW[i] = 0;
		L_a_BW[i] = INF;
		beta[i] = 0;
	}
	for (int j = 0; j < DEVICES; j++) {
		I_edge[j].clear();
		I_cloud[j].clear();
		f_cloud[j] = 0;
		energy_comp_cloud[j] = 0;
		eta[j] = 0;
		device[j].p = 0;
		device[j].task.r_BW = 0;
		device[j].task.inputData = 0;
		device[j].task.outputData = 0;
		for (int k = 0; k < RESOURCE; k++) {
			device[j].task.r[k] = 0;
		}
		for (int i = 0; i < ACCESS; i++) {
			device[j].SNR[i] = 0;
			device[j].h[i] = 0;   //信道增益
			device[j].v[i] = 0;

			energy_MT[i][j] = 0;  
			energy_edgeserver2[i][j] = 0;   
			energy_tran2[i][j] = 0;
		}
	}
	cloud.constrain_BW = 0;
	cloud.z_BW = 0;
	for (int k = 0; k < RESOURCE; k++) {
		cloud.constrain_r[k] = 0;
		cloud.z_r[k] = 0;
		U_c[k] = 0;
		L_c[k] = INF;
		gama[k] = 0;

	}
	U_c_BW = 0;
	L_c_BW = INF;
	deta = 0;

	start_init = 0;
	end_init = 0;
}

// 从文件中读取实验数据
void read_data(string fdevice, string faccess, string fSNR) {
	//根据file中数据初始化数据
	string data;
	ifstream infdevice;
	ifstream infaccess;
	ifstream infSNR;

	infdevice.open(fdevice);
	//设备数据，内存大小，输入输出数据大小都是MB
	for (int j = 0; j < DEVICES; j++)
	{
		getline(infdevice, data);//按行读入 
		istringstream istr1(data); //从string对象data中读取空格、换行符分割的字符
		for (int k = 0; k < RESOURCE; k++)
		{
			istr1 >> device[j].task.r[k];    //内存和disk的单位是MB
		}

		istr1 >> device[j].task.inputData;
		device[j].task.inputData *= 1;  //MB
		istr1 >> device[j].task.outputData;
		device[j].task.outputData *= 1;   //MB  
		istr1 >> device[j].task.r_BW;
		device[j].task.r_BW *= 1;            //MHz


		cout << device[j].task.r[0] << '\t' << device[j].task.r[1] << '\t' << device[j].task.r[2] << '\t' << device[j].task.inputData << endl;

		istr1 >> f_cloud[j];
		f_cloud[j] *= 1;  //GHz

		istr1 >> device[j].p;  //0.2-0.5w

		getline(infdevice, data);
		istringstream istr_f(data);
		for (int i = 0; i < ACCESS; i++) {
			istr_f >> f_edge[i][j];
			f_edge[i][j] *= 1;    //GHz
		}
	}
	infdevice.close();

	infaccess.open(faccess);
	//无线接入点
	for (int i = 0; i < ACCESS; i++)
	{
		getline(infaccess, data);
		istringstream istr2(data);
		for (int k = 0; k < RESOURCE; k++)
		{
			istr2 >> access[i].constrain_r[k];
			access[i].constrain_r[k] *= 1;
		}
		access[i].constrain_r[0] *= 1;  //Gcycles 
		access[i].constrain_r[1] *= 1;    //MB
		access[i].constrain_r[2] *= 1;    //MB
		istr2 >> access[i].constrain_BW;
		access[i].constrain_BW *= 1;    //MHz
		istr2 >> access[i].p;                //w
		istr2 >> access[i].v;                //MBps

		cout << access[i].constrain_r[0] << '\t' << access[i].constrain_r[1] << '\t' << access[i].constrain_r[2] << '\t' << access[i].constrain_BW << '\t' << access[i].p << '\t' << access[i].v << endl;
	}


	infaccess.close();

	infSNR.open(fSNR);
	//信噪比
	for (int j = 0; j < DEVICES; j++)
	{
		getline(infSNR, data);
		istringstream istr3(data);
		for (int i = 0; i < ACCESS; i++) {
			istr3 >> device[j].SNR[i];
		}
	}
	infSNR.close();

	//云端
	
	cloud.constrain_r[0] = 6000;  
	cloud.constrain_r[1] = 1024 * 12; 
	cloud.constrain_r[2] = 1024 * 10;	
	cloud.constrain_BW = 3000; 

}

// 计算对偶变量
void calculate_price(char p, int i, int k) {
	
	double index;
	switch (p) {
	case 'a': {
		//计算alpha
		index = access[i].z_r[k] / access[i].constrain_r[k];
		alpha[i][k] = L_a[i][k] * pow(U_a[i][k] / L_a[i][k], index);
	}; break;
	case 'g': {
		//计算gama
		index = cloud.z_r[k] / cloud.constrain_r[k];
		gama[k] = L_c[k] * pow(U_c[k] / L_c[k], index);
	}; break;
	case 'b': {
		//计算beta
		index = access[i].z_BW / access[i].constrain_BW;
		beta[i] = L_a_BW[i] * pow(U_a_BW[i] / L_a_BW[i], index);
	}; break;
	case 'd': {
		//计算deta
		index = (double)cloud.z_BW / cloud.constrain_BW;
		//cout << "c.z_BW = " << c.z_BW << " c.constrain_BW = " << c.constrain_BW << " index = " << index;
		deta = L_c_BW * pow(U_c_BW / L_c_BW, index);
		//cout << "deta:" << deta << endl;
	}; break;
	default:break;
	}
}

// 初始化U、L的值
void init_UL() {
	double temp;
	for (int i = 0; i < ACCESS; i++) {
		for (int k = 0; k < RESOURCE; k++) {
			for (int j = 1; j < DEVICES; j++) {
				//cout << temp << " ";
				if (k == 0)
					temp = Energy_edge[i][j] / (device[j].task.r[k] * 1024 * 1024 * 1024);  
				else
					temp = Energy_edge[i][j] / (device[j].task.r[k] * 1024 * 1024 * 8);
				if (U_a[i][k] < temp) U_a[i][k] = temp;
				if (L_a[i][k] > temp) L_a[i][k] = temp;
			}
		}
		for (int j = 1; j < DEVICES; j++) {
			temp = Energy_edge[i][j] / ((double)device[j].task.r_BW * 1024 * 1024);  
			if (U_a_BW[i] < temp) U_a_BW[i] = temp;
			if (L_a_BW[i] > temp) L_a_BW[i] = temp;
		}


		for (int k = 0; k < RESOURCE; k++) {
			for (int j = 0; j < DEVICES; j++) {
				if (k == 0)
					temp = Energy_cloud[i][j] / (device[j].task.r[k] * 1024 * 1024 * 1024);  
				else
					temp = Energy_cloud[i][j] / (device[j].task.r[k] * 1024 * 1024 * 8);
				if (U_c[k] < temp) U_c[k] = temp;
				if (L_c[k] > temp) L_c[k] = temp;
			}

		}
		for (int j = 0; j < DEVICES; j++) {
			temp = Energy_cloud[i][j] / ((double)device[j].task.r_BW * 1024 * 1024); 
			if (U_c_BW < temp) U_c_BW = temp;
			if (L_c_BW > temp) L_c_BW = temp;
		}
	}
}

// 判断边缘资源是否溢出
bool is_full_cpubw(int i, int j) {    //两种资源
	for (int k = 0; k < RESOURCE; k++) {
		if (access[i].z_r[k] + device[j].task.r[k] > access[i].constrain_r[k])
			return false;
	}
	if (access[i].z_BW + device[j].task.r_BW > access[i].constrain_BW)
		return false;
	return true;
}

bool is_full_bw(int i, int j) {    //BW资源
	if (access[i].z_BW + device[j].task.r_BW > access[i].constrain_BW)
		return false;
	return true;
}

// 判断云资源是否溢出
bool is_full(int j) {
	for (int k = 0; k < RESOURCE; k++) {
		if (cloud.z_r[k] + device[j].task.r[k] > cloud.constrain_r[k])
			return false;
	}
	if (cloud.z_BW + device[j].task.r_BW > cloud.constrain_BW)
		return false;
	return true;
}

void redefinition()
{
	for (int i = 0; i < ACCESS; i++)
	{
		for (int k = 0; k < RESOURCE; k++) {
			if (k == 0)//su
				U_a[i][k] = U_a[i][k] * 1024 * 1024 * 1024;//su
			else//su
				U_a[i][k] = U_a[i][k] * 1024 * 1024 * 8;//su
			L_a[i][k] = U_a[i][k] / UdL;
		}
		U_a_BW[i] = U_a_BW[i] * 1024 * 1024;//su
		L_a_BW[i] = U_a_BW[i] / UdL;
	}

	for (int k = 0; k < RESOURCE; k++) {
		if (k == 0)
			U_c[k] = U_c[k] * 1024 * 1024 * 1024;//su
		else
			U_c[k] = U_c[k] * 1024 * 1024 * 8;//su
		L_c[k] = U_c[k] / UdL;
	}
	U_c_BW = U_c_BW * 1024 * 1024;//su
	L_c_BW = U_c_BW / UdL;

}

// 计算能耗
void cal_energy() {

	for (int i = 0; i < ACCESS; i++) {
		for (int j = 0; j < DEVICES; j++) {

			device[j].v[i] = device[j].task.r_BW * (log(1 + device[j].SNR[i]) / log(2));//单位：Mbps            
			//迁移到边缘计算
			energy_tran_access[i][j] = (device[j].p / 1e-3) * (device[j].task.inputData * 8 / device[j].v[i] / 3600) * 3.6 * 1e6;   //单位：j
			energy_comp_edge[i][j] = theta1 * device[j].task.r[0] * 1024 * 1024 * 1024 * f_edge[i][j] * 1024 * 1024 * 1024 * f_edge[i][j] * 1024 * 1024 * 1024;   //单位：j   
			energy_tran_cloud1[i][j] = access[i].p / 1e-3 * (device[j].task.outputData * 8 / access[i].v / 3600) * 3.6 * 1e6 + device[j].task.outputData * ENERGY_INTENSITY; //单位：j
			//cout << i << "," << j << endl;
			//cout << access[i].p << '\t' << device[j].task.outputData << '\t' << access[i].v << endl;
			//cout << energy_tran_cloud1[i][j] << endl << endl;
			Energy_edge[i][j] = energy_tran_access[i][j] + energy_comp_edge[i][j] + energy_tran_cloud1[i][j];

			energy_edgeserver1[i][j] = energy_comp_edge[i][j] + access[i].p / 1e-3 * (device[j].task.outputData / access[i].v / 3600) * 3.6 * 1e6;    //迁移到边的边缘服务器能耗 
			energy_tran1[i][j] = energy_tran_access[i][j] + energy_tran_cloud1[i][j]; //迁移到边的传输能耗  

			//迁移到云端计算
			energy_tran_cloud2[i] = access[i].p / 1e-3 * (device[j].task.inputData * 8 / access[i].v / 3600) * 3.6 * 1e6 + device[j].task.inputData * ENERGY_INTENSITY;  //单位：j
			energy_comp_cloud[j] = theta2 * device[j].task.r[0] * 1024 * 1024 * 1024 * f_cloud[j] * 1024 * 1024 * 1024 * f_cloud[j] * 1024 * 1024 * 1024;  //单位：j
			Energy_cloud[i][j] = energy_tran_access[i][j] + energy_tran_cloud2[i] + energy_comp_cloud[j];

			energy_MT[i][j] = energy_tran_access[i][j];  //迁移到云的用户终端能耗 
			energy_edgeserver2[i][j] = access[i].p / 1e-3 * (device[j].task.inputData * 8 / access[i].v / 3600) * 3.6 * 1e6;  //迁移到云的边缘服务器能耗 
			energy_tran2[i][j] = energy_tran_access[i][j] + energy_tran_cloud2[i];//迁移到云的传输能耗 
		}
	}

}

// 更新sigma
void update_sigma(int i)
{
	double fenzi = 0;
	double fenmu = 0;
	for (int j = 0; j < DEVICES; j++)
	{
		fenzi += x[i][j] * device[j].task.outputData;
		fenmu += x[i][j] * device[j].task.inputData;
	}
	access[i].sigma = fenzi / fenmu;
}

// 算法主体
void EB_CMA(int j) {
	double rjk = 0;

	for (int i = 0; i < ACCESS; i++) {
		if (is_full_bw(i, j))
			I_cloud[j].push_back(i);
		if (is_full_cpubw(i, j))
			I_edge[j].push_back(i);
	}

	for (int i = 0; i < I_cloud[j].size(); i++) {
		rjk = 0;
		for (int k = 0; k < RESOURCE; k++) {
			rjk += device[j].task.r[k] * gama[k];
		}
		eta_cloud[I_cloud[j][i]][j] = Energy_cloud[I_cloud[j][i]][j] + (device[j].task.r_BW * beta[I_cloud[j][i]] + rjk + device[j].task.r_BW * deta);
	}

	rjk = 0;

	for (int i = 0; i < I_edge[j].size(); i++) {
		rjk = 0;
		for (int k = 0; k < RESOURCE; k++) {
			rjk += device[j].task.r[k] * alpha[I_edge[j][i]][k];
		}
		eta_edge[I_edge[j][i]][j] = Energy_edge[I_edge[j][i]][j] + (rjk + device[j].task.r_BW * beta[I_edge[j][i]] + access[I_edge[j][i]].sigma * device[j].task.r_BW * deta); //access[I_edge[j][i]].sigma 
	}

	eta[j] = INF;
	int i_cloud = -1;  //记录最小的eta_cloud[i][j]是经过哪个i中继求得的

	int i0 = -2;		//标记，-1表示迁移到云端计算，-2表示不迁移。若i0 >= 0, 则迁移到无线接入点，i0是无线接入点的编号


	for (int i = 0; i < I_cloud[j].size(); i++) {
		if (eta_cloud[I_cloud[j][i]][j] < eta[j])
		{
			eta[j] = eta_cloud[I_cloud[j][i]][j];
			i_cloud = I_cloud[j][i];
			i0 = -1;
		}

	}

	// min{eta_edge}
	for (int i = 0; i < I_edge[j].size(); i++) {
		if (eta_edge[I_edge[j][i]][j] < eta[j])
		{
			eta[j] = eta_edge[I_edge[j][i]][j];
			i0 = I_edge[j][i];
		}
	}

	if (i0 >= 0) {
		// 分配到边，边缘服务器下标为i0
		x[i0][j] = 1;
		update_sigma(i0);
		for (int k = 0; k < RESOURCE; k++) {
			access[i0].z_r[k] += device[j].task.r[k];
			calculate_price('a', i0, k);
		}

		access[i0].z_BW += device[j].task.r_BW;
		cloud.z_BW += access[i0].sigma * device[j].task.r_BW;

		calculate_price('b', i0, 0);
		calculate_price('d', 0, 0);
	}
	else if (i0 == -1) {
		// 分配到云
		y[i_cloud][j] = 1;     

		access[i_cloud].z_BW += device[j].task.r_BW;
		calculate_price('b', i_cloud, 0);

		for (int k = 0; k < RESOURCE; k++) {
			cloud.z_r[k] += device[j].task.r[k];
			calculate_price('g', 0, k);
		}
		cloud.z_BW += device[j].task.r_BW;
		calculate_price('d', 0, 0);
	}

	I_cloud.clear();
	I_edge.clear();
}

void print_result(Result* cplex_r, Result* cplexcloud_r, Result* cplexedge_r, Result* cma_r, Result* ff_r, Result* e_r)  //在控制台显示
{
	/* 需要打印的结果数据包括：
	*  1.用户分配到边计算的比例
	*  2.总能量消耗
	*  3.云/边各种资源的使用率
	*  4.时间消耗
	* */

	cout << "------------------------------结果------------------------------" << endl;
	cout << "\t\tCplex\t\tCplexcloud\tCplexedge\tPrimalDual\tFirstFit\tAverage\n";
	// 1.用户分配到边计算的比例
	if (isPrint.p_edge_uti == 1) {
		cout << "edge_uti\t" << setprecision(5) << cplex_r->edge_uti << "\t\t" << setprecision(5) << cplexcloud_r->edge_uti << "\t\t" << setprecision(5) << cplexedge_r->edge_uti << "\t\t" << setprecision(5) << cma_r->edge_uti << "\t\t";
		cout << setprecision(5) << ff_r->edge_uti << "\t\t" << setprecision(5) << e_r->edge_uti << "\n";
	}
	// 2.总能量消耗
	if (isPrint.p_result_energy == 1) {
		cout << "energy\t\t" << setprecision(12) << cplex_r->energy << "\t" << setprecision(12) << cplexcloud_r->energy << "\t" << setprecision(12) << cplexedge_r->energy << "\t" << setprecision(12) << cma_r->energy << "\t";
		cout << setprecision(12) << ff_r->energy << "\t" << setprecision(12) << e_r->energy << "\n";
	}
	// 3.云/边各种资源的使用率
	if (isPrint.p_uti == 1) {
		for (int k = 0; k < RESOURCE; k++) {
			cout << "r" << k + 1 << "_uti\t\t" << setprecision(5) << cplex_r->uti[k] << "\t\t" << setprecision(5) << cplexcloud_r->uti[k] << "\t\t" << setprecision(5) << cplexedge_r->uti[k] << "\t\t" << setprecision(5) << cma_r->uti[k] << "\t\t";
			cout << setprecision(5) << ff_r->uti[k] << "\t\t" << setprecision(5) << e_r->uti[k] << "\n";
		}
		cout << "BW_uti\t\t" << setprecision(5) << cplex_r->uti_BW << "\t\t" << setprecision(5) << cplexcloud_r->uti_BW << "\t\t" << setprecision(5) << cplexedge_r->uti_BW << "\t\t" << setprecision(5) << cma_r->uti_BW << "\t\t";
		cout << setprecision(5) << ff_r->uti_BW << "\t\t" << setprecision(5) << e_r->uti_BW << "\n";
	}
	if (isPrint.p_uti_c == 1) {
		for (int k = 0; k < RESOURCE; k++) {
			cout << "r" << k + 1 << "_uti_c\t" << setprecision(5) << cplex_r->uti_c[k] << "\t" << setprecision(5) << cplexcloud_r->uti_c[k] << "\t" << setprecision(5) << cplexedge_r->uti_c[k] << "\t\t" << setprecision(5) << cma_r->uti_c[k] << "\t";
			cout << setprecision(5) << ff_r->uti_c[k] << "\t" << setprecision(5) << e_r->uti_c[k] << "\n";
		}
		cout << setprecision(5) << "BW_uti_c\t" << cplex_r->uti_c_BW << "\t\t" << setprecision(5) << cplexcloud_r->uti_c_BW << "\t\t" << setprecision(5) << cplexedge_r->uti_c_BW << "\t\t" << setprecision(5) << cma_r->uti_c_BW << "\t\t";
		cout << setprecision(5) << ff_r->uti_c_BW << "\t\t" << setprecision(5) << e_r->uti_c_BW << "\n";
	}
	// 4.时间
	if (isPrint.p_t == 1) {
		cout << "time\t\t" << setprecision(12) << cplex_r->t << "\t\t" << setprecision(12) << cplexcloud_r->t << "\t\t" << setprecision(12) << cplexedge_r->t << "\t\t" << setprecision(12) << cma_r->t << "\t\t";
		cout << setprecision(12) << ff_r->t << "\t\t" << setprecision(12) << e_r->t << "\n";
	}
	//5.用户终端、边缘服务器能耗及占比  
	if (isPrint.p_energy_MT == 1) {
		cout << "en_MT\t\t" << setprecision(12) << cplex_r->energy_MT << "\t" << setprecision(12) << cplexcloud_r->energy_MT << "\t" << setprecision(12) << cplexedge_r->energy_MT << "\t" << setprecision(12) << cma_r->energy_MT << "\t";
		cout << setprecision(12) << ff_r->energy_MT << "\t" << setprecision(12) << e_r->energy_MT << "\n";
	}
	if (isPrint.p_energy_edgeservercom == 1) {
		cout << "en_edgecom\t" << setprecision(12) << cplex_r->energy_edgeservercom << "\t" << setprecision(12) << cplexcloud_r->energy_edgeservercom << "\t" << setprecision(12) << cplexedge_r->energy_edgeservercom << "\t" << setprecision(12) << cma_r->energy_edgeservercom << "\t";
		cout << setprecision(12) << ff_r->energy_edgeservercom << "\t" << setprecision(12) << e_r->energy_edgeservercom << "\n";
	}
	if (isPrint.p_energy_edgeserver == 1) {
		cout << "en_edge\t\t" << setprecision(12) << cplex_r->energy_edgeserver << "\t" << setprecision(12) << cplexcloud_r->energy_edgeserver << "\t" << setprecision(12) << cplexedge_r->energy_edgeserver << "\t" << setprecision(12) << cma_r->energy_edgeserver << "\t";
		cout << setprecision(12) << ff_r->energy_edgeserver << "\t" << setprecision(12) << e_r->energy_edgeserver << "\n";
	}
	if (isPrint.p_energy_tran == 1) {
		cout << "en_tran\t\t" << setprecision(12) << cplex_r->energy_tran << "\t" << setprecision(12) << cplexcloud_r->energy_tran << "\t" << setprecision(12) << cplexedge_r->energy_tran << "\t" << setprecision(12) << cma_r->energy_tran << "\t";
		cout << setprecision(12) << ff_r->energy_tran << "\t" << setprecision(12) << e_r->energy_tran << "\n";
	}

}

void print_csv(Result* cplex_r, Result* cplexcloud_r, Result* cplexedge_r, Result* cma_r, Result* ff_r, Result* e_r)  //输出到CSV文件
{
	// 需要打印的结果数据包括：
	//  1.用户分配到边计算的比例
	// 2.总能量消耗
	// 3.云/边各种资源的使用率
	// 4.时间消耗

	ofstream output;
   string root = "CSV/numofMT/sigma0.25/";
	//string root = "CSV/tasktype/sigma0.25/";
	//string root = "CSV/edgeMTratio/";
	//string root = "CSV/sigma/";
	//string root = "CSV/UL/sigma0.8/";
	string ename = "energy.csv";
	string urname;
	string ubname = "uti_bw.csv"; //边带宽利用率
	string ubcname = "uti_bwcloud.csv"; //云带宽利用率
	string uename = "uti_edge.csv";  //任务分配到边计算的比例
	string tname = "time.csv";

	string eMTname = "energyMT.csv";  
	string eedgesercomname = "energyedgesercom.csv";   
	string eedgesername = "energyedgeser.csv";   
	string etranname = "energytran.csv";
	

	// 1.用户分配到边计算的比例
	if (isPrint.p_edge_uti == 1) {
		output.open(root + uename, ios::app);
		output << setprecision(5) << cplex_r->edge_uti << ',' << setprecision(5) << cplexcloud_r->edge_uti << ',' << setprecision(5) << cplexedge_r->edge_uti << ',' << setprecision(5) << cma_r->edge_uti << ',';
		output << setprecision(5) << ff_r->edge_uti << ',' << setprecision(5) << e_r->edge_uti << endl;
		output.close();
	}
	// 2.总能量消耗
	if (isPrint.p_result_energy == 1) {
		output.open(root + ename, ios::app);
		output << setprecision(12) << cplex_r->energy << ',' << setprecision(12) << cplexcloud_r->energy << ',' << setprecision(12) << cplexedge_r->energy << ',' << setprecision(12) << cma_r->energy << ',';
		output << setprecision(12) << ff_r->energy << ',' << setprecision(12) << e_r->energy << "\n";
		output.close();
	}
	// 3.云/边各种资源的使用率
	if (isPrint.p_uti == 1) {
		for (int k = 0; k < RESOURCE; k++) {
			urname = "ur" + to_string(k) + "edge.csv";
			output.open(root + urname, ios::app);
			output << setprecision(5) << cplex_r->uti[k] << ',' << setprecision(5) << cplexcloud_r->uti[k] << ',' << setprecision(5) << cplexedge_r->uti[k] << ',' << setprecision(5) << cma_r->uti[k] << ',';
			output << setprecision(5) << ff_r->uti[k] << ',' << setprecision(5) << e_r->uti[k] << "\n";
			output.close();
		}
		output.open(root + ubname, ios::app);
		output << setprecision(5) << cplex_r->uti_BW << ',' << setprecision(5) << cplexcloud_r->uti_BW << ',' << setprecision(5) << cplexedge_r->uti_BW << ',' << setprecision(5) << cma_r->uti_BW << ',';
		output << setprecision(5) << ff_r->uti_BW << ',' << setprecision(5) << e_r->uti_BW << "\n";
		output.close();

		for (int k = 0; k < RESOURCE; k++) {
			urname = "ur" + to_string(k) + "cloud.csv";
			output.open(root + urname, ios::app);
			output << setprecision(5) << cplex_r->uti_c[k] << ',' << setprecision(5) << cplexcloud_r->uti_c[k] << ',' << setprecision(5) << cplexedge_r->uti_c[k] << ',' << setprecision(5) << cma_r->uti_c[k] << ',';
			output << setprecision(5) << ff_r->uti_c[k] << ',' << setprecision(5) << e_r->uti_c[k] << "\n";
			output.close();
		}

		output.open(root + ubcname, ios::app);
		output << setprecision(5) << cplex_r->uti_c_BW << ',' << setprecision(5) << cplexcloud_r->uti_c_BW << ',' << setprecision(5) << cplexedge_r->uti_c_BW << ',' << setprecision(5) << cma_r->uti_c_BW << ',';
		output << setprecision(5) << ff_r->uti_c_BW << ',' << setprecision(5) << e_r->uti_c_BW << "\n";
		output.close();
	}
	// 4.时间
	if (isPrint.p_t == 1) {
		output.open(root + tname, ios::app);
		output << setprecision(12) << cplex_r->t << ',' << setprecision(12) << cplexcloud_r->t << ',' << setprecision(12) << cplexedge_r->t << ',' << setprecision(12) << cma_r->t << ',';
		output << setprecision(12) << ff_r->t << ',' << setprecision(12) << e_r->t << "\n";
		output.close();
	}
	//5.用户终端、边缘服务器能耗及占比  su
	if (isPrint.p_energy_MT == 1) {
		output.open(root + eMTname, ios::app);
		output << setprecision(12) << cplex_r->energy_MT << ',' << setprecision(12) << cplexcloud_r->energy_MT << ',' << setprecision(12) << cplexedge_r->energy_MT << ',' << setprecision(12) << cma_r->energy_MT << ',';
		output << setprecision(12) << ff_r->energy_MT << ',' << setprecision(12) << e_r->energy_MT << endl;
		output.close();
	}
	if (isPrint.p_energy_edgeservercom == 1) {
		output.open(root + eedgesercomname, ios::app);
		output << setprecision(12) << cplex_r->energy_edgeservercom << ',' << setprecision(12) << cplexcloud_r->energy_edgeservercom << ',' << setprecision(12) << cplexedge_r->energy_edgeservercom << ',' << setprecision(12) << cma_r->energy_edgeservercom << ',';
		output << setprecision(12) << ff_r->energy_edgeservercom << ',' << setprecision(12) << e_r->energy_edgeservercom << endl;
		output.close();
	}
	if (isPrint.p_energy_edgeserver == 1) {
		output.open(root + eedgesername, ios::app);
		output << setprecision(12) << cplex_r->energy_edgeserver << ',' << setprecision(12) << cplexcloud_r->energy_edgeserver << ',' << setprecision(12) << cplexedge_r->energy_edgeserver << ',' << setprecision(12) << cma_r->energy_edgeserver << ',';
		output << setprecision(12) << ff_r->energy_edgeserver << ',' << setprecision(12) << e_r->energy_edgeserver << endl;
		output.close();
	}
	if (isPrint.p_energy_tran == 1) {
		output.open(root + etranname, ios::app);
		output << setprecision(12) << cplex_r->energy_tran << ',' << setprecision(12) << cplexcloud_r->energy_tran << ',' << setprecision(12) << cplexedge_r->energy_tran << ',' << setprecision(12) << cma_r->energy_tran << ',';
		output << setprecision(12) << ff_r->energy_tran << ',' << setprecision(12) << e_r->energy_tran << "\n";
		output.close();
	}
}

// 每组实验开始，初始化参数（设置初值、从文件读取数据）
void init(string fdevice, string faccess, string fSNR) {
	//初始化
	for (int i = 0; i < ACCESS; i++)
	{
		for (int k = 0; k < RESOURCE; k++)
			L_a[i][k] = INF;
		L_a_BW[i] = INF;
	}
	for (int k = 0; k < RESOURCE; k++)
		L_c[k] = INF;
	L_c_BW = INF;

	start_init = clock();
	//根据test_case所指的文件，初始化设备与无线接入点的初始信息
	read_data(fdevice, faccess, fSNR);
	//根据已知信息计算能耗
	cal_energy();
	end_init = clock();
	//根据已知信息，计算U与L的值
	init_UL();

	//根据计算好的UL值，初始化边际价格
	for (int i = 0; i < ACCESS; i++) {
		calculate_price('b', i, 0);
		for (int k = 0; k < RESOURCE; k++) {
			calculate_price('a', i, k);
		}
	}
	for (int k = 0; k < RESOURCE; k++) {
		calculate_price('d', 0, k);
	}
	calculate_price('e', 0, 0);
}

// 打印初始参数
void print_data()
{
	/*isPrint.p_xy = 0;
	isPrint.p_initial_condition = 0;
	isPrint.p_now_energy[0] = 0;
	isPrint.p_now_energy[1] = 0;
	isPrint.p_UL = 0;
	isPrint.p_dual = 0;
	isPrint.p_miu = 0;
	isPrint.p_I = 0;*/
	int row = 0;
	if (ACCESS > 15)
		row = 15;
	else
		row = ACCESS;
	if (isPrint.p_initial_condition == 1) {
		cout << "\n------------------------------cloud data------------------------------" << endl;
		cout << "tc1\tc2\tc3\tc_BW\tz1\tz2\tz3\tz_BW\n";
		cout << setprecision(5) << cloud.constrain_r[0] << "\t" << cloud.constrain_r[1] << "\t" << cloud.constrain_r[2] << "\t" << cloud.constrain_BW << "\t";
		cout << setprecision(5) << cloud.z_r[0] << "\t" << cloud.z_r[1] << "\t" << cloud.z_r[2] << "\t" << cloud.z_BW << '\n';

		cout << "\n------------------------------devices data------------------------------" << endl;
		cout << "\tr1\tr2\tr3\tBW\tinput\toutput\tp\tSNR\tcpu_f\n";
		for (int j = 0; j < DEVICES; j++)
		{
			cout << "d" << j + 1 << "\t";
			for (int k = 0; k < RESOURCE; k++)
			{
				cout << setprecision(5) << device[j].task.r[k] << '\t';
			}

			cout << setprecision(5) << device[j].task.r_BW << '\t';
			cout << setprecision(5) << device[j].task.inputData << '\t';
			cout << setprecision(5) << device[j].task.outputData << '\t';
			cout << device[j].p << '\t';
			cout << device[j].SNR[0] << '\t';
		}

		cout << "\n------------------------------accesses data------------------------------" << endl;
		cout << "\t" << "tr1\t\tz1\t\tr2\t\tz2\t\tr3\t\tz3\t\tBW\tz_BW\tp\tv\tsigma\n";
		for (int i = 0; i < ACCESS; i++)
		{
			cout << "a" << i + 1 << "\t";
			for (int k = 0; k < RESOURCE; k++)
			{
				cout << setprecision(6) << access[i].constrain_r[k] << "\t\t";
				cout << setprecision(6) << access[i].z_r[k] << "\t\t";
			}
			cout << access[i].constrain_BW << '\t';
			cout << access[i].z_BW << '\t';
			cout << access[i].p << '\t';
			cout << access[i].v << '\t';
			cout << setprecision(12) << access[i].sigma << endl;
		}
	}


	if (isPrint.p_now_energy[0] == 1) {
		cout << "\n------------------------------devices to access energy------------------------------" << endl;
		cout << '\t';
		for (int i = 0; i < row; i++)
		{
			cout << "a" << i + 1 << "\t\t";
		}
		cout << '\n';
		for (int j = 0; j < 30; j++)
		{
			cout << "d" << j + 1 << '\t';
			for (int i = 0; i < row; i++)
			{
				cout << setprecision(6) << energy_tran_access[i][j] << '\t';
			}
			cout << endl;
		}

		cout << "\n------------------------------energy of comp at access------------------------------" << endl;
		cout << '\t';
		for (int i = 0; i < row; i++)
		{
			cout << "a" << i + 1 << "\t";
		}
		cout << '\n';
		for (int j = 0; j < 30; j++)
		{
			cout << "d" << j + 1 << '\t';
			for (int i = 0; i < row; i++)
			{
				cout << setprecision(6) << energy_comp_edge[i][j] << '\t';
			}
			cout << endl;
		}

		cout << "\n------------------------------energy of output tran to cloud------------------------------" << endl;
		cout << '\t';
		for (int i = 0; i < row; i++)
		{
			cout << "a" << i + 1 << "\t";
		}
		cout << '\n';
		for (int j = 0; j < 30; j++)
		{
			cout << "d" << j + 1 << '\t';
			for (int i = 0; i < row; i++)
			{
				cout << setprecision(6) << energy_tran_cloud1[i][j] << '\t';
			}
			cout << endl;
		}

		cout << "\n------------------------------total energy of comp at access------------------------------" << endl;
		cout << '\t';
		for (int i = 0; i < row; i++)
		{
			cout << "a" << i + 1 << "\t";
		}
		cout << '\n';
		for (int j = 0; j < 30; j++)
		{
			cout << "d" << j + 1 << '\t';
			for (int i = 0; i < row; i++)
			{
				cout << setprecision(6) << Energy_edge[i][j] << '\t';
			}
			cout << endl;
		}

		cout << "\n------------------------------total energy of comp at cloud------------------------------" << endl;
		for (int j = 0; j < 30; j++)
		{
			cout << "d" << j + 1 << " ";
		}
		cout << '\n';
		for (int j = 0; j < 30; j++)
		{
			cout << Energy_cloud[j] << '\t';
		}
		cout << endl;
	}


	if (isPrint.p_UL == 1) {
		cout << endl << "U_a(ik):" << endl;
		for (int i = 0; i < ACCESS; i++)
		{
			for (int k = 0; k < RESOURCE; k++)
				cout << U_a[i][k] << '\t';
			cout << endl;
		}

		cout << endl << "L_a(ik):" << endl;
		for (int i = 0; i < ACCESS; i++)
		{
			for (int k = 0; k < RESOURCE; k++)
				cout << L_a[i][k] << '\t';
			cout << endl;
		}

		cout << endl << "U_a_BW(i):" << endl;
		for (int i = 0; i < ACCESS; i++)
		{
			cout << U_a_BW[i] << '\t';
			if (i + 1 % 10 == 0)
				cout << '\n';
		}

		cout << endl;
		cout << endl << "L_a_BW(i):" << endl;
		for (int i = 0; i < ACCESS; i++) {
			cout << L_a_BW[i] << '\t';
			if (i + 1 % 10 == 0)
				cout << '\n';
		}

		cout << endl;
		cout << endl << "U_c(k):" << endl;
		for (int k = 0; k < RESOURCE; k++)
			cout << U_c[k] << '\t';

		cout << endl;
		cout << endl << "L_c(k):" << endl;
		for (int k = 0; k < RESOURCE; k++)
			cout << L_c[k] << '\t';

		cout << endl;
		cout << endl << "U_c_BW:" << endl;
		cout << U_c_BW << '\t';

		cout << endl;
		cout << endl << "L_c_BW:" << endl;
		cout << L_c_BW << '\t';

		cout << endl << endl;

		/// <summary>
		/// 
		/// </summary>
		double* a_ = new double[RESOURCE] {0};                  
		double a_bw = 0;
		double temp = 0;
		cout << endl << "U_a(ik)/L_a(ik):" << endl;
		for (int i = 0; i < ACCESS; i++)
		{
			for (int k = 0; k < RESOURCE; k++) {
				temp = U_a[i][k] / L_a[i][k];
				a_[k] += temp;
				cout << temp << '\t';
			}
			cout << endl;
		}
		for (int k = 0; k < RESOURCE; k++) {
			ULa[k] += a_[k] / ACCESS;    //边上的U/L,取平均
			cout << a_[k] / ACCESS << '\t';
		}
		cout << endl;

		cout << endl << "U_a_BW(i)/L_a_BW(i):" << endl;
		for (int i = 0; i < ACCESS; i++)
		{
			temp = U_a_BW[i] / L_a_BW[i];
			a_bw += temp;
			cout << temp << '\t';
			if (i + 1 % 10 == 0)
				cout << '\n';
		}
		ULBWa += a_bw / ACCESS;
		cout << '\n' << a_bw / ACCESS << '\n';

		cout << endl;
		cout << endl << "U_c(k)/L_c(k):" << endl;
		for (int k = 0; k < RESOURCE; k++) {
			temp = U_c[k] / L_c[k];
			ULc[k] += temp;
			cout << temp << '\t';
		}

		cout << endl;
		cout << endl << "U_c_BW/L_c_BW:" << endl;
		temp = U_c_BW / L_c_BW;
		ULBWc += temp;
		cout << temp << '\t';
		cout << endl << endl;
	}

	if (isPrint.p_dual == 1) {
		cout << "arpha:" << endl;
		for (int i = 0; i < ACCESS; i++) {
			for (int k = 0; k < RESOURCE; k++) {
				cout << alpha[i][k] << " ";
			}
			cout << endl;
		}
		cout << endl << endl;

		cout << "gama:" << endl;
		for (int k = 0; k < RESOURCE; k++) {
			cout << gama[k] << " ";
		}
		cout << endl << endl;

		cout << "beta:" << endl;
		for (int i = 0; i < ACCESS; i++) {
			cout << beta[i] << "\t";
			if (i + 1 % 10 == 0)
				cout << '\n';
		}
		cout << endl << endl;

		cout << "deta:" << deta << endl << endl;
	}

	if (isPrint.p_miu == 1) {
		cout << "eta:" << endl;
		for (int j = 0; j < 30; j++) {
			cout << eta[j] << "\t";
			if (j + 1 % 10 == 0)
				cout << endl;
		}
		cout << endl << endl;

		cout << "eta_cloud:" << endl;
		for (int j = 0; j < 30; j++) {
			cout << eta_cloud[j] << "\t";
			if (j + 1 % 10 == 0)
				cout << endl;
		}
		cout << endl << endl;

		cout << "eta_edge:" << endl;
		for (int j = 0; j < 30; j++) {
			for (int i = 0; i < 10; i++) {
				cout << eta_edge[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << endl;
	}


}


// 优先分配边
void first_fit_alg(Result* result)
{
	//按照first fit方法来决策
	//需要用到的全局变量有：
	//1.每个设备的请求资源数 d[j].t.r[k] d[j].t.r_BW[k] 
	//2.每个无线接入点的资源限制 a[i].constrain_r a[i].constrain_BW
	//3.每种策略所消耗的能量 Energy_cloud[j] Energy_edge[i][j]
	//需要声明的新变量有
	//1.每个无线接入点当前的被分配的资源  r[i][k] r_BW[i]
	//2.总的能量消耗Energy_ff
	int(*x_ff)[DEVICES] = new int[ACCESS][DEVICES]{ 0 };
	int(*y_ff)[DEVICES] = new int[ACCESS][DEVICES]{ 0 };
	double(*z)[RESOURCE] = new double[ACCESS][RESOURCE]{ 0 };
	double* z_BW = new double[ACCESS] { 0 };
	double* z_c = new double[RESOURCE] { 0 };  //云上已分配
	double z_BW_c = 0;       //云上已分配
	//标记变量count，标记被分配的无线接入点编号，若当前分配成功，则值为非负
	int count = -1;
	for (int j = 0; j < DEVICES; j++) {
		for (int i = 0; i < ACCESS; i++) {
			count = -1;
			if (device[j].task.r_BW + z_BW[i] > access[i].constrain_BW)
				continue;
			for (int k = 0; k < RESOURCE; k++) {
				if (device[j].task.r[k] + z[i][k] > access[i].constrain_r[k])
				{
					count = -2;
					break;
				}
			}

			if (count == -2)
				continue;
			else
			{
				//
				count = i;
				for (int k = 0; k < RESOURCE; k++)
					z[i][k] += device[j].task.r[k];
				z_BW[i] += device[j].task.r_BW;
				z_BW_c += device[j].task.r_BW * 0.25;
				x_ff[i][j] = 1;
				break;
			}
		}
		if (count < 0) {
			/*if (device[j].task.r_BW + r_BW_c > cloud.constrain_BW)
				continue;
			for (int k = 0; k < RESOURCE; k++) {
				if (device[j].task.r[k] + r_c[k] > cloud.constrain_r[k])
				{
					count = -1;
					break;
				}
			}
			if (count == -1)
				continue;*/
		//所有无线接入点都无资源分配，计算迁移至云，任意选一个access作中继
			for (int i = 0; i < ACCESS; i++)
				if (device[j].task.r_BW + z_BW[i] <= access[i].constrain_BW)
				{
					y_ff[i][j] = 1;
					z_BW[i] += device[j].task.r_BW;
					break;
				}

			for (int k = 0; k < RESOURCE; k++)
				z_c[k] += device[j].task.r[k];
			z_BW_c += device[j].task.r_BW;
		}
	}

	cal_result(result, x_ff, y_ff, isPrint.p_xy);


		delete[] z_c;
		delete[] z_BW;
		delete[] z;
		delete[] y_ff;
		delete[] x_ff;
	}
}

// 平均分配
void equal_alg(Result* result)
{
	//按照均匀迁移方法来决策
	//需要用到的全局变量有：
	//1.每个设备的请求资源数 d[j].t.r[k] d[j].t.r_BW[k] 
	//2.每个无线接入点的资源限制 a[i].constrain_r a[i].constrain_BW
	//3.每种策略所消耗的能量 Energy_cloud[j] Energy_edge[i][j]
	//需要声明的新变量有
	//1.每个无线接入点当前的被分配的资源  r[i][k] r_BW[i]
	//2.总的能量消耗Energy_ff
	int(*x_e)[DEVICES] = new int[ACCESS][DEVICES]{ 0 };  
	int(*y_e)[DEVICES] = new int[ACCESS][DEVICES]{ 0 };
	double(*z)[RESOURCE] = new double[ACCESS][RESOURCE]{ 0 };
	double* z_BW = new double[ACCESS] { 0 };
	double* z_c = new double[RESOURCE] { 0 };
	double z_BW_c = 0;
	//标记变量count，标记被分配的无线接入点编号，若当前分配成功，则值为非负
	int count = -1;
	for (int j = 0; j < DEVICES; j++) {
		if (j % 2 == 0) {
			for (int i = 0; i < ACCESS; i++) {
				count = -1;
				if (device[j].task.r_BW + z_BW[i] > access[i].constrain_BW)
					continue;
				for (int k = 0; k < RESOURCE; k++) {
					if (device[j].task.r[k] + z[i][k] > access[i].constrain_r[k])
					{
						count = -2;
						break;
					}
				}

				if (count == -2)
					continue;
				else
				{
					count = i;
					for (int k = 0; k < RESOURCE; k++)
						z[i][k] += device[j].task.r[k];
					z_BW[i] += device[j].task.r_BW;
					z_BW_c += device[j].task.r_BW * 0.25;
					x_e[i][j] = 1;
					break;
				}
			}
			if (count < 0) {
				/*if (device[j].task.r_BW + r_BW_c > cloud.constrain_BW)
					continue;
				for (int k = 0; k < RESOURCE; k++) {
					if (device[j].task.r[k] + r_c[k] > cloud.constrain_r[k])
					{
						count = -1;
						break;
					}
				}
				if (count == -1)
					continue;*/
				//所有无线接入点都无资源分配，计算迁移至云
				for (int i = 0; i < ACCESS; i++)
					if (device[j].task.r_BW + z_BW[i] <= access[i].constrain_BW)
					{
						y_e[i][j] = 1;
						z_BW[i] += device[j].task.r_BW;
						break;
					}

				for (int k = 0; k < RESOURCE; k++)
					z_c[k] += device[j].task.r[k];
				z_BW_c += device[j].task.r_BW;
			}
		}
		else
		{

			for (int i = 0; i < ACCESS; i++)
				if (device[j].task.r_BW + z_BW[i] <= access[i].constrain_BW)
				{
					y_e[i][j] = 1;
					z_BW[i] += device[j].task.r_BW;
					break;
				}

			for (int k = 0; k < RESOURCE; k++)
				z_c[k] += device[j].task.r[k];
			z_BW_c += device[j].task.r_BW;
		}
	}

	cal_result(result, x_e, y_e, isPrint.p_xy);


	delete[] z_c;
	delete[] z_BW;
	delete[] z;
	delete[] y_e;
	delete[] x_e;

}



void rukou() {
	Result cplex_r;
	Result cplexcloud_r;
	Result cplexedge_r;
	Result cma_r;
	Result ff_r;
	Result e_r;

	cplex_r = *new Result;
	cplexcloud_r = *new Result;
	cplexedge_r = *new Result;
	cma_r = *new Result;
	ff_r = *new Result;
	e_r = *new Result;

	//时间变量
	clock_t start_1, start_2, start_3, start_4, start_5, start_6, end_1, end_2, end_3, end_4, end_5, end_6;
	int d_num = 1;		//每种情况的实验组数
	int index = 0;

	//UdL =5;   // U/L=1.2,1.5,2,exp(1),10,20,50,80,100

	// 打印参数
	isPrint.p_xy = 0;
	isPrint.p_initial_condition = 0;
	isPrint.p_now_energy[0] = 0;
	isPrint.p_now_energy[1] = 0;
	isPrint.p_UL = 1;
	isPrint.p_dual = 0;
	isPrint.p_miu = 0;
	//isPrint.p_I = 0;

	isPrint.p_result_energy = 1;
	isPrint.p_energy_MT = 1; 
	isPrint.p_energy_edgeservercom = 1; 
	isPrint.p_energy_edgeserver = 1; 
	isPrint.p_energy_tran = 1;

	isPrint.p_edge_uti = 1;
	isPrint.p_uti = 1;
	isPrint.p_uti_c = 1;
	isPrint.p_t = 1;


	//设置数据源路径
	string root = "F:/实验202208/DATA/numofMT/sigma0.25/500/";
	//string root = "F:/实验202208/DATA/tasktype/sigma0.25/100/";
	//string root="C:/Users/SQ/Desktop/实验202208/DATA/sigma/1/";
	//string root="C:/Users/SQ/Desktop/实验202208/DATA/edgeMTratio/20/";
	//string root = "C:/Users/SQ/Desktop/实验202208/DATA/UL/sigma1/";
	for (int n = 0; n < d_num; n++) {

		// 每次实验都重置一次全局变量(重置为0)
		reset_all();


		/// <summary>
		/// 构造实验数据文件路径
		/// </summary>
		string fdevice = root + "data_" + to_string(DEVICES) + "devices" + to_string(n + index) + ".txt";	// 终端数据
		string faccess = root + "data_" + to_string(int(DEVICES * 0.04)) + "access" + to_string(n + index) + ".txt";		// 边缘接入点数据
		string fSNR = root + "data_" + to_string(int(DEVICES * 0.04)) + "access" + to_string(DEVICES) + "devices" + to_string(n + index) + ".txt";		// 信噪比

		/// <summary>
		/// 根据数据文件初始化所有已知量


		start_1 = clock();
		init(fdevice, faccess, fSNR);


	   //redefinition();


		//本文方法
		// 一个一个地分配
		for (int j = 0; j < DEVICES; j++)
			EB_CMA(j);
		cal_result(&cma_r, x, y, isPrint.p_xy);
		end_1 = clock();
		cma_r.t += (double)(end_1)-(double)(start_1);
		print_data();

		//first_fit
		start_2 = clock();
		first_fit_alg(&ff_r);
		end_2 = clock();
		ff_r.t += ((double)(end_init)-(double)(start_init)) + ((double)(end_2)-(double)(start_2));

		//equal
		start_3 = clock();
		equal_alg(&e_r);
		end_3 = clock();
		e_r.t += ((double)(end_init)-(double)(start_init)) + ((double)(end_3)-(double)(start_3));
		//cout << "etime:" << e_time <<  endl;

		//cplex
		start_4 = clock();
		////if (n < d_num) {
		cplex(&cplex_r);
		end_4 = clock();
		cplex_r.t += ((double)(end_init)-(double)(start_init)) + ((double)(end_4)-(double)(start_4));

		//cplex_cloud
		start_5 = clock();
		cplex_cloud(&cplexcloud_r);
		end_5 = clock();
		cplexcloud_r.t += ((double)(end_init)-(double)(start_init)) + ((double)(end_5)-(double)(start_5));
		//cout << cplexcloud_r.t << "-----------------[" << endl;

		//cplex_edge
		start_6 = clock();
		cplex_edge(&cplexedge_r);
		end_6 = clock();
		cplexedge_r.t += ((double)(end_init)-(double)(start_init)) + ((double)(end_6)-(double)(start_6));
		//cout << cplexedge_r.t << "-----------------[" << endl;

	}


	divide_num(&cplex_r, d_num);
	divide_num(&cplexcloud_r, d_num);
	divide_num(&cplexedge_r, d_num);
	divide_num(&cma_r, d_num);
	divide_num(&ff_r, d_num);
	divide_num(&e_r, d_num);

	cout << "--------------UL---------------\n";
	cout << "a:\t";
	for (int k = 0; k < RESOURCE; k++)
		cout << ULa[k] / d_num << '\t';
	cout << ULBWa / d_num << '\n';
	cout << "c:\t";
	for (int k = 0; k < RESOURCE; k++)
		cout << ULc[k] / d_num << '\t';
	cout << ULBWc / d_num << '\n' << '\n';


	print_result(&cplex_r, &cplexcloud_r, &cplexedge_r, &cma_r, &ff_r, &e_r);
	print_csv(&cplex_r, &cplexcloud_r, &cplexedge_r, &cma_r, &ff_r, &e_r);
}

void runtime()
{
	Result cplex_r;
	Result cplexcloud_r;
	Result cplexedge_r;
	Result cma_r;
	Result ff_r;
	Result e_r;

	cplex_r = *new Result;
	cplexcloud_r = *new Result;
	cplexedge_r = *new Result;
	cma_r = *new Result;
	ff_r = *new Result;
	e_r = *new Result;

	//时间变量
	clock_t start_1, start_2, start_3, start_4, start_5, start_6, end_1, end_2, end_3, end_4, end_5, end_6;
	//int d_num = 100;		//每种情况的实验组数
	int index = 0;

	//UdL =5;   // U/L=1.2,1.5,2,exp(1),10,20,50,80,100

	// 打印参数
	isPrint.p_xy = 1;
	isPrint.p_initial_condition = 0;
	isPrint.p_now_energy[0] = 0;
	isPrint.p_now_energy[1] = 0;
	isPrint.p_UL = 1;
	isPrint.p_dual = 0;
	isPrint.p_miu = 0;
	//isPrint.p_I = 0;

	isPrint.p_result_energy = 1;
	isPrint.p_energy_MT = 1; 
	isPrint.p_energy_edgeservercom = 1; 
	isPrint.p_energy_edgeserver = 1; 
	isPrint.p_energy_tran = 1;

	isPrint.p_edge_uti = 1;
	isPrint.p_uti = 1;
	isPrint.p_uti_c = 1;
	isPrint.p_t = 1;


	//设置数据源路径
	string root = "F:/实验202208/DATA/numofMT/sigma0.25/500/";
	//string root = "C:/Users/SQ/Desktop/实验202208/DATA/tasktype/100/";
	//string root="C:/Users/SQ/Desktop/实验202208/DATA/sigma/1/";
	//string root="C:/Users/SQ/Desktop/实验202208/DATA/edgeMTratio/20/";
	//string root = "C:/Users/SQ/Desktop/实验202208/DATA/UL/sigma1/";
	//for (int n = 0; n < d_num; n++) {


	// 每次实验都重置一次全局变量(重置为0)
	reset_all();


	/// <summary>
	/// 构造实验数据文件路径
	/// </summary>
	string fdevice = root + "data_" + to_string(DEVICES) + "devices" + to_string(0 + index) + ".txt";	// 终端数据
	string faccess = root + "data_" + to_string(int(DEVICES * 0.04)) + "access" + to_string(0 + index) + ".txt";		// 边缘接入点数据
	string fSNR = root + "data_" + to_string(int(DEVICES * 0.04)) + "access" + to_string(DEVICES) + "devices" + to_string(0 + index) + ".txt";		// 信噪比

	/// <summary>
	/// 根据数据文件初始化所有已知量
	/// </summary>

	start_1 = clock();
	init(fdevice, faccess, fSNR);

   //redefinition();

	//cplex
	start_4 = clock();
	cplex(&cplex_r);
	end_4 = clock();
	cplex_r.t += ((double)(end_init)-(double)(start_init)) + ((double)(end_4)-(double)(start_4));


	//cplex_cloud
	start_5 = clock();
	cplex_cloud(&cplexcloud_r);
	end_5 = clock();
	cplexcloud_r.t += ((double)(end_init)-(double)(start_init)) + ((double)(end_5)-(double)(start_5));

	//cplex_edge
	start_6 = clock();
	cplex_edge(&cplexedge_r);
	end_6 = clock();
	cplexedge_r.t += ((double)(end_init)-(double)(start_init)) + ((double)(end_6)-(double)(start_6));


	//本文方法
	// 一个一个地分配
	for (int j = 0; j < DEVICES; j++) {
		cma_r.energy = 0;
		EB_CMA(j);
		cal_result(&cma_r, x, y, isPrint.p_xy);
		print_result(&cplex_r, &cplexcloud_r, &cplexedge_r, &cma_r, &ff_r, &e_r);
		print_csv(&cplex_r, &cplexcloud_r, &cplexedge_r, &cma_r, &ff_r, &e_r);
	}

	//end_1 = clock();
	//cma_r.t += (double)(end_1)-(double)(start_1);
	//print_data();

	//first_fit
	/*start_2 = clock();
	first_fit_alg(&ff_r);
	end_2 = clock();
	ff_r.t += ((double)(end_init)-(double)(start_init)) + ((double)(end_2)-(double)(start_2));

	//equal
	start_3 = clock();
	equal_alg(&e_r);
	end_3 = clock();
	e_r.t += ((double)(end_init)-(double)(start_init)) + ((double)(end_3)-(double)(start_3));
	//cout << "etime:" << e_time <<  endl;*/

	//}


	/*divide_num(&cplex_r, d_num);
	divide_num(&cplexcloud_r, d_num);
	divide_num(&cplexedge_r, d_num);
	divide_num(&cma_r, d_num);
	divide_num(&ff_r, d_num);
	divide_num(&e_r, d_num);

	cout << "--------------UL---------------\n";
	cout << "a:\t";
	for (int k = 0; k < RESOURCE; k++)
		cout << ULa[k] / d_num << '\t';
	cout << ULBWa / d_num << '\n';
	cout << "c:\t";
	for (int k = 0; k < RESOURCE; k++)
		cout << ULc[k] / d_num << '\t';
	cout << ULBWc / d_num << '\n' << '\n';*/


	//print_result(&cplex_r, &cplexcloud_r, &cplexedge_r, &cma_r, &ff_r, &e_r);
	//print_csv(&cplex_r, &cplexcloud_r, &cplexedge_r, &cma_r, &ff_r, &e_r);
}

int main()
{
	rukou();
	//runtime();

}
