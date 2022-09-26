#include <vector>
#include <limits.h>
#define INF LONG_MAX  //Maximum value for a variable of type long:2147483647
#define RESOURCE 3

#define DEVICES 100
#define ACCESS 4


//能量强度，单位 焦耳/MB
#define ENERGY_INTENSITY (0.06* 3.6 * 1e6) / 1024  //0.06kWh/GB，单位：j/MB,1千瓦时=3600000焦  


// 是否打印变量的控制符
static struct IsPrint {
	int p_xy = 0;
	int p_initial_condition = 0;
	int p_now_energy[2] = { 0 };
	int p_UL = 0;
	int p_dual = 0;
	int p_miu = 0;
	int p_I = 0;

	int p_result_energy = 0;
	int p_energy_MT = 0;  
	int p_energy_edgeservercom = 0; 
	int p_energy_edgeserver = 0; 
	int p_energy_tran = 0; //总传输能耗 

	int p_edge_uti = 0;
	int p_uti = 0;
	int p_uti_c = 0;
	int p_t = 0;

} isPrint;

// 保存结果的结构体
struct Result {

	int x[ACCESS][DEVICES] = { 0 }; //原始变量
	int y[ACCESS][DEVICES] = { 0 };  //原始变量 
	double energy = 0; // 能耗

	//用户终端、边缘服务器能耗
	double energy_MT = 0;  
	double energy_edgeservercom = 0;  
	double energy_edgeserver = 0;  
	double energy_tran = 0;//总传输能耗  

	double edge_uti = 0; // 任务分配到边的比例
	double uti[RESOURCE] = { 0 }; // 各计算资源利用率
	double uti_BW = { 0 };  //带宽资源利用率
	double uti_c[RESOURCE] = { 0 };  //云上各计算资源利用率
	double uti_c_BW = 0;

	double t = 0; // 时间消耗

	 //用户终端、边缘服务器能耗在总能耗中的占比
	double energy_MT_uti = 0;   
	double energy_edgeservercom_uti = 0;   
	double energy_edgeserver_uti = 0;  
	double energy_tran_uti = 0;//总传输能耗占比  
};

struct Resource {
	int cpu;
	int memory;
	int disk;
};

struct Task {
	double r[RESOURCE] = { 0 };
	int r_BW = 0;
	double inputData;
	double outputData;
};

//无线终端
struct Devices {
	int SNR[ACCESS];	//设备j到无线接入点i的信噪比
	double h[ACCESS];		//设备j到无线接入点i的信道增益
	double v[ACCESS];		//设备j到无线接入点i的信道最大传输速率
	double p;		//设备j的信号发射功率
	Task task;
};

//无线接入点
struct Access {
	double constrain_r[RESOURCE] = { 0.0 };   //资源容量
	int constrain_BW = 0;  //带宽资源容量
	double sigma;  //边缘服务器上输出数据量与输入数据量之比
	double v;    //互联网传输速率
	double p;   //无线接入点的发射功率
	double z_r[RESOURCE] = { 0 };  //无线接入点已经被分配掉的计算资源
	int z_BW = 0;  //无线接入点已经被分配掉的带宽资源
};

struct Cloud {
	double constrain_r[RESOURCE] = { 0.0 };
	int constrain_BW = 0;
	double z_r[RESOURCE] = { 0 };
	int z_BW = 0;
};

//以上结构体都是数据类型定义，下面才是具体的数据对象

Devices device[DEVICES];

Access access[ACCESS];

Cloud cloud;

//theta i，与无线接入点i硬件结构相关的常数
double theta1 = 1e-26;
double theta2 = 1e-26;
//信道噪声N[i][j]
double N[ACCESS][DEVICES] = { 0.0 };

int x[ACCESS][DEVICES] = { 0 };
int y[ACCESS][DEVICES] = { 0 };  



double f_edge[ACCESS][DEVICES] = { 0 };		//i为j分配的算力
double f_cloud[DEVICES] = { 0 };				//云为j分配的算力

double time_tran[ACCESS][DEVICES] = { 0.0 };  //终端到无线接入点的传输时间

double energy_tran_access[ACCESS][DEVICES] = { 0.0 };			//迁移到边，终端传输到边缘的能耗
double energy_comp_edge[ACCESS][DEVICES] = { 0.0 };				//在边的计算能耗
double energy_comp_cloud[DEVICES] = { 0.0 };					//在云的计算能耗
double energy_tran_cloud1[ACCESS][DEVICES] = { 0.0 };			//迁移到边，输出数据从边到云备份的传输能耗
double energy_tran_cloud2[ACCESS] = { 0.0 };					//迁移到云，从边到云的传输能耗

double Energy_edge[ACCESS][DEVICES] = { 0.0 };					//迁移到边的总能耗
double Energy_cloud[ACCESS][DEVICES] = { 0.0 };					//迁移到云的总能耗 

double energy_edgeserver1[ACCESS][DEVICES] = { 0.0 }; //迁移到边的边缘服务器总能耗 
double energy_tran1[ACCESS][DEVICES] = { 0.0 }; //迁移到边的总传输能耗 

double energy_MT[ACCESS][DEVICES] = { 0.0 }; //迁移到云的用户终端能耗  
double energy_edgeserver2[ACCESS][DEVICES] = { 0.0 };  //迁移到云的边缘服务器能耗 
double energy_tran2[ACCESS][DEVICES] = { 0.0 };//迁移到云的总传输能耗 

/// <summary>
/// 
/// </summary>
double U_a[ACCESS][RESOURCE] = { 0.0 };    //边上各计算资源的最大价值
double U_a_BW[ACCESS] = { 0.0 };          //边上带宽资源的最大价值
double L_a[ACCESS][RESOURCE] = { INF };  //边上各计算资源的最小价值
double L_a_BW[ACCESS] = { INF };    //边上带宽资源的最小价值

double U_c[RESOURCE] = { 0.0 };    //云上各计算资源的最大价值
double U_c_BW = { 0.0 };
double L_c[RESOURCE] = { INF };
double L_c_BW = { INF };


/// <summary>
/// 对偶变量
/// </summary>
double alpha[ACCESS][RESOURCE] = { 0.0 };
double gama[RESOURCE] = { 0.0 };
double beta[ACCESS] = { 0.0 };
double deta = 0.0;
double eta[DEVICES] = { 0.0 };  //应该是epsilon
double eta_cloud[ACCESS][DEVICES] = { 0.0 };  
double eta_edge[ACCESS][DEVICES] = { 0.0 };

double ULa[RESOURCE] = { 0 };   //U/L
double ULBWa = 0;
double ULc[RESOURCE] = { 0 };
double ULBWc = 0;

double UdL = 0;  // U/L

clock_t start_init = 0;
clock_t end_init = 0;

std::vector<std::vector<int>> I_edge(DEVICES);  //相当于二维动态数组，其中第一个维度为DEVICES
std::vector<std::vector<int>> I_cloud(DEVICES);  


void cal_result(Result* result, const int(*x0)[DEVICES], const int(*y0)[DEVICES], int p_xy)  
{
	// 初始化当前组的计算结果
	double now_energy = 0;

	double now_energy_MT = 0;   
	double now_energy_edgeservercom = 0;  
	double now_energy_edgeserver = 0;  
	double now_energy_tran = 0;

	double now_uti[RESOURCE] = { 0 };  //当前边上的计算资源利用率
	double now_uti_BW = { 0 };
	double now_uti_c[RESOURCE] = { 0 };   //当前云上的计算资源利用率
	double now_uti_c_BW = 0;
	double z_r[ACCESS][RESOURCE] = { 0 };   //z表示已分配资源量
	double z_BW[ACCESS] = { 0 };
	double z_c_r[RESOURCE] = { 0 };
	double z_c_BW = 0;


	// 计算分配到边的能耗及资源使用情况
	for (int i = 0; i < ACCESS; i++) {
		for (int j = 0; j < DEVICES; j++) {
			now_energy += x0[i][j] * Energy_edge[i][j];
			now_energy_MT += x0[i][j] * energy_tran_access[i][j];  
			now_energy_edgeservercom += x0[i][j] * energy_comp_edge[i][j];  
			now_energy_edgeserver += x0[i][j] * energy_edgeserver1[i][j];   
			now_energy_tran += x0[i][j] * energy_tran1[i][j];

			for (int k = 0; k < RESOURCE; k++) {
				z_r[i][k] += (double)x0[i][j] * device[j].task.r[k];
			}
			z_BW[i] += ((double)x0[i][j] * device[j].task.r_BW + (double)y0[i][j] * device[j].task.r_BW);   
		}
	}
	for (int i = 0; i < ACCESS; i++) {
		for (int k = 0; k < RESOURCE; k++) {
			now_uti[k] += (double)z_r[i][k] / access[i].constrain_r[k];
		}
		now_uti_BW += (double)z_BW[i] / access[i].constrain_BW;
	}
	for (int k = 0; k < RESOURCE; k++) {
		result->uti[k] += now_uti[k] / ACCESS;   //平均利用率
	}
	result->uti_BW += now_uti_BW / ACCESS;  //平均利用率

	// 计算分配到云的能耗及资源使用情况
	for (int i = 0; i < ACCESS; i++) {
		for (int j = 0; j < DEVICES; j++) {
			now_energy += y0[i][j] * Energy_cloud[i][j];
			now_energy_MT += y0[i][j] * energy_MT[i][j];  
			now_energy_edgeserver += y0[i][j] * energy_edgeserver2[i][j];  
			now_energy_tran += y0[i][j] * energy_tran2[i][j];

			for (int k = 0; k < RESOURCE; k++) {
				z_c_r[k] += y0[i][j] * device[j].task.r[k];
			}
			z_c_BW += ((double)y0[i][j] * device[j].task.r_BW + (double)x0[i][j] * device[j].task.r_BW * access[i].sigma);   
		}
	}
	result->energy += now_energy;

	result->energy_MT += now_energy_MT;  
	result->energy_edgeservercom += now_energy_edgeservercom;  
	result->energy_edgeserver += now_energy_edgeserver;  
	result->energy_tran += now_energy_tran; 

	//计算云端的资源利用率
	for (int k = 0; k < RESOURCE; k++) {
		now_uti_c[k] = z_c_r[k] / cloud.constrain_r[k];
		result->uti_c[k] += now_uti_c[k];
	}
	now_uti_c_BW = z_c_BW / cloud.constrain_BW;
	result->uti_c_BW += now_uti_c_BW;

	double* count = new double;
	*count = 0;
	for (int i = 0; i < ACCESS; i++) {
		for (int j = 0; j < DEVICES; j++) {
			if (x0[i][j] == 1)
				*count = *count + 1;   //分配到边执行的任务数
		}
	}
	*count /= (double)DEVICES;   //分配到边执行的任务比例
	result->edge_uti += *count;

	if (isPrint.p_xy == 1) {
		for (int i = 0; i < ACCESS; i++) {
			for (int j = 0; j < DEVICES; j++) {
				cout << x0[i][j] << " ";
			}
			cout << '\n';
		}
		cout << '\n';
		for (int i = 0; i < ACCESS; i++) {
			for (int j = 0; j < DEVICES; j++) {
				cout << y0[i][j] << " ";
			}
			cout << '\n';
			cout << '\n';
		}
	}

}

void divide_num(Result* result, int num) {    //执行多次取平均值
	result->edge_uti /= (double)num;
	result->energy /= (double)num;

	result->energy_MT /= (double)num;  
	result->energy_edgeservercom /= (double)num;  
	result->energy_edgeserver /= (double)num;  
	result->energy_tran /= (double)num;

	result->t /= (double)num;
	for (int k = 0; k < RESOURCE; k++) {
		result->uti[k] /= (double)num;
		result->uti_c[k] /= (double)num;
	}

	result->uti_BW /= (double)num;
	result->uti_c_BW /= (double)num;
}