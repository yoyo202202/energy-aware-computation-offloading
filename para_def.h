#include <vector>
#include <limits.h>
#define INF LONG_MAX  //Maximum value for a variable of type long:2147483647
#define RESOURCE 3

#define DEVICES 100
#define ACCESS 4


//����ǿ�ȣ���λ ����/MB
#define ENERGY_INTENSITY (0.06* 3.6 * 1e6) / 1024  //0.06kWh/GB����λ��j/MB,1ǧ��ʱ=3600000��  


// �Ƿ��ӡ�����Ŀ��Ʒ�
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
	int p_energy_tran = 0; //�ܴ����ܺ� 

	int p_edge_uti = 0;
	int p_uti = 0;
	int p_uti_c = 0;
	int p_t = 0;

} isPrint;

// �������Ľṹ��
struct Result {

	int x[ACCESS][DEVICES] = { 0 }; //ԭʼ����
	int y[ACCESS][DEVICES] = { 0 };  //ԭʼ���� 
	double energy = 0; // �ܺ�

	//�û��նˡ���Ե�������ܺ�
	double energy_MT = 0;  
	double energy_edgeservercom = 0;  
	double energy_edgeserver = 0;  
	double energy_tran = 0;//�ܴ����ܺ�  

	double edge_uti = 0; // ������䵽�ߵı���
	double uti[RESOURCE] = { 0 }; // ��������Դ������
	double uti_BW = { 0 };  //������Դ������
	double uti_c[RESOURCE] = { 0 };  //���ϸ�������Դ������
	double uti_c_BW = 0;

	double t = 0; // ʱ������

	 //�û��նˡ���Ե�������ܺ������ܺ��е�ռ��
	double energy_MT_uti = 0;   
	double energy_edgeservercom_uti = 0;   
	double energy_edgeserver_uti = 0;  
	double energy_tran_uti = 0;//�ܴ����ܺ�ռ��  
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

//�����ն�
struct Devices {
	int SNR[ACCESS];	//�豸j�����߽����i�������
	double h[ACCESS];		//�豸j�����߽����i���ŵ�����
	double v[ACCESS];		//�豸j�����߽����i���ŵ����������
	double p;		//�豸j���źŷ��书��
	Task task;
};

//���߽����
struct Access {
	double constrain_r[RESOURCE] = { 0.0 };   //��Դ����
	int constrain_BW = 0;  //������Դ����
	double sigma;  //��Ե�����������������������������֮��
	double v;    //��������������
	double p;   //���߽����ķ��书��
	double z_r[RESOURCE] = { 0 };  //���߽�����Ѿ���������ļ�����Դ
	int z_BW = 0;  //���߽�����Ѿ���������Ĵ�����Դ
};

struct Cloud {
	double constrain_r[RESOURCE] = { 0.0 };
	int constrain_BW = 0;
	double z_r[RESOURCE] = { 0 };
	int z_BW = 0;
};

//���Ͻṹ�嶼���������Ͷ��壬������Ǿ�������ݶ���

Devices device[DEVICES];

Access access[ACCESS];

Cloud cloud;

//theta i�������߽����iӲ���ṹ��صĳ���
double theta1 = 1e-26;
double theta2 = 1e-26;
//�ŵ�����N[i][j]
double N[ACCESS][DEVICES] = { 0.0 };

int x[ACCESS][DEVICES] = { 0 };
int y[ACCESS][DEVICES] = { 0 };  



double f_edge[ACCESS][DEVICES] = { 0 };		//iΪj���������
double f_cloud[DEVICES] = { 0 };				//��Ϊj���������

double time_tran[ACCESS][DEVICES] = { 0.0 };  //�ն˵����߽����Ĵ���ʱ��

double energy_tran_access[ACCESS][DEVICES] = { 0.0 };			//Ǩ�Ƶ��ߣ��ն˴��䵽��Ե���ܺ�
double energy_comp_edge[ACCESS][DEVICES] = { 0.0 };				//�ڱߵļ����ܺ�
double energy_comp_cloud[DEVICES] = { 0.0 };					//���Ƶļ����ܺ�
double energy_tran_cloud1[ACCESS][DEVICES] = { 0.0 };			//Ǩ�Ƶ��ߣ�������ݴӱߵ��Ʊ��ݵĴ����ܺ�
double energy_tran_cloud2[ACCESS] = { 0.0 };					//Ǩ�Ƶ��ƣ��ӱߵ��ƵĴ����ܺ�

double Energy_edge[ACCESS][DEVICES] = { 0.0 };					//Ǩ�Ƶ��ߵ����ܺ�
double Energy_cloud[ACCESS][DEVICES] = { 0.0 };					//Ǩ�Ƶ��Ƶ����ܺ� 

double energy_edgeserver1[ACCESS][DEVICES] = { 0.0 }; //Ǩ�Ƶ��ߵı�Ե���������ܺ� 
double energy_tran1[ACCESS][DEVICES] = { 0.0 }; //Ǩ�Ƶ��ߵ��ܴ����ܺ� 

double energy_MT[ACCESS][DEVICES] = { 0.0 }; //Ǩ�Ƶ��Ƶ��û��ն��ܺ�  
double energy_edgeserver2[ACCESS][DEVICES] = { 0.0 };  //Ǩ�Ƶ��Ƶı�Ե�������ܺ� 
double energy_tran2[ACCESS][DEVICES] = { 0.0 };//Ǩ�Ƶ��Ƶ��ܴ����ܺ� 

/// <summary>
/// 
/// </summary>
double U_a[ACCESS][RESOURCE] = { 0.0 };    //���ϸ�������Դ������ֵ
double U_a_BW[ACCESS] = { 0.0 };          //���ϴ�����Դ������ֵ
double L_a[ACCESS][RESOURCE] = { INF };  //���ϸ�������Դ����С��ֵ
double L_a_BW[ACCESS] = { INF };    //���ϴ�����Դ����С��ֵ

double U_c[RESOURCE] = { 0.0 };    //���ϸ�������Դ������ֵ
double U_c_BW = { 0.0 };
double L_c[RESOURCE] = { INF };
double L_c_BW = { INF };


/// <summary>
/// ��ż����
/// </summary>
double alpha[ACCESS][RESOURCE] = { 0.0 };
double gama[RESOURCE] = { 0.0 };
double beta[ACCESS] = { 0.0 };
double deta = 0.0;
double eta[DEVICES] = { 0.0 };  //Ӧ����epsilon
double eta_cloud[ACCESS][DEVICES] = { 0.0 };  
double eta_edge[ACCESS][DEVICES] = { 0.0 };

double ULa[RESOURCE] = { 0 };   //U/L
double ULBWa = 0;
double ULc[RESOURCE] = { 0 };
double ULBWc = 0;

double UdL = 0;  // U/L

clock_t start_init = 0;
clock_t end_init = 0;

std::vector<std::vector<int>> I_edge(DEVICES);  //�൱�ڶ�ά��̬���飬���е�һ��ά��ΪDEVICES
std::vector<std::vector<int>> I_cloud(DEVICES);  


void cal_result(Result* result, const int(*x0)[DEVICES], const int(*y0)[DEVICES], int p_xy)  
{
	// ��ʼ����ǰ��ļ�����
	double now_energy = 0;

	double now_energy_MT = 0;   
	double now_energy_edgeservercom = 0;  
	double now_energy_edgeserver = 0;  
	double now_energy_tran = 0;

	double now_uti[RESOURCE] = { 0 };  //��ǰ���ϵļ�����Դ������
	double now_uti_BW = { 0 };
	double now_uti_c[RESOURCE] = { 0 };   //��ǰ���ϵļ�����Դ������
	double now_uti_c_BW = 0;
	double z_r[ACCESS][RESOURCE] = { 0 };   //z��ʾ�ѷ�����Դ��
	double z_BW[ACCESS] = { 0 };
	double z_c_r[RESOURCE] = { 0 };
	double z_c_BW = 0;


	// ������䵽�ߵ��ܺļ���Դʹ�����
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
		result->uti[k] += now_uti[k] / ACCESS;   //ƽ��������
	}
	result->uti_BW += now_uti_BW / ACCESS;  //ƽ��������

	// ������䵽�Ƶ��ܺļ���Դʹ�����
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

	//�����ƶ˵���Դ������
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
				*count = *count + 1;   //���䵽��ִ�е�������
		}
	}
	*count /= (double)DEVICES;   //���䵽��ִ�е��������
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

void divide_num(Result* result, int num) {    //ִ�ж��ȡƽ��ֵ
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