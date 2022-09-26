#pragma once
#include <ilcplex/ilocplex.h>
#include <stdio.h>
#include <iostream>
using namespace std;

#include "para_def.h"

typedef IloArray <IloNumVarArray> IloNumVarArray2;

ILOSTLBEGIN


//double Energy_cloud_cplex[DEVICES] = { 0 };
//double Energy_edge_cplex[ACCESS][DEVICES] = { 0 };
//double C_edge[ACCESS][RESOURCE] = { 0 };
//double C_edge_BW[ACCESS] = { 0 };
//double C_cloud[RESOURCE] = { 0 };
//double C_cloud_BW = 0;
//double r[DEVICES][RESOURCE] = { 0 };
//double r_BW[DEVICES] = { 0 };
//double output_data[DEVICES] = { 0 };
//double input_data[DEVICES] = { 0 };

//最小化云、边能耗
int xc[ACCESS][DEVICES] = { 0 };
int yc[ACCESS][DEVICES] = { 0 };

//最小化云能耗
int xc_cloud[ACCESS][DEVICES] = { 0 };
int yc_cloud[ACCESS][DEVICES] = { 0 };

//最小化边能耗
int xc_edge[ACCESS][DEVICES] = { 0 };
int yc_edge[ACCESS][DEVICES] = { 0 };



//最小化云、边能耗
void cplex(Result* result)
{

	IloEnv env;
	double solution_value = 0;
	try {

		IloModel model(env);

		IloNumVarArray2 x(env, ACCESS);
		IloNumVarArray2 y(env, ACCESS);
		for (IloInt i = 0; i < ACCESS; i++)
		{
			x[i] = IloNumVarArray(env, DEVICES);
			y[i] = IloNumVarArray(env, DEVICES);
			for (IloInt j = 0; j < DEVICES; j++) {
				x[i][j] = IloNumVar(env, 0, 1, ILOINT);
				model.add(x[i][j]);

				y[i][j] = IloNumVar(env, 0, 1, ILOINT);
				model.add(y[i][j]);
			}
		}


		//约束8a
		IloExpr sum_edge_r(env);
		for (IloInt i = 0; i < ACCESS; i++)
		{
			for (IloInt k = 0; k < RESOURCE; k++)
			{
				sum_edge_r.clear();
				for (IloInt j = 0; j < DEVICES; j++)
					sum_edge_r += device[j].task.r[k] * x[i][j];
				model.add(sum_edge_r <= access[i].constrain_r[k]);
			}
		}


		//约束8b
		IloExpr sum_edge_BW(env);
		for (IloInt i = 0; i < ACCESS; i++)
		{
			sum_edge_BW.clear();
			for (IloInt j = 0; j < DEVICES; j++)
			{
				sum_edge_BW += device[j].task.r_BW * (x[i][j] + y[i][j]);
			}
			model.add(sum_edge_BW <= access[i].constrain_BW);
		}

		//约束8c
		IloExpr sum_cloud_r(env);
		for (IloInt k = 0; k < RESOURCE; k++)
		{
			sum_cloud_r.clear();
			for (IloInt i = 0; i < ACCESS; i++)
			{
				for (IloInt j = 0; j < DEVICES; j++)
				{
					sum_cloud_r += device[j].task.r[k] * y[i][j];
				}
			}
			model.add(sum_cloud_r <= cloud.constrain_r[k]);
		}

		//约束8d
		IloExpr sum_cloud_BW(env);
		for (IloInt i = 0; i < ACCESS; i++)
		{
			for (IloInt j = 0; j < DEVICES; j++) {
				sum_cloud_BW += device[j].task.r_BW * y[i][j];
			}
		}
		IloExpr temp(env);
		for (IloInt i = 0; i < ACCESS; i++)
		{
			temp.clear();
			for (IloInt j = 0; j < DEVICES; j++) {
				temp += device[j].task.r_BW * x[i][j];
			}
			try {
				sum_cloud_BW += access[i].sigma * temp;		
			}
			catch (exception e) {
				sum_cloud_BW += 0;
			}
		}
		model.add(sum_cloud_BW <= cloud.constrain_BW);

		//约束8e
		IloExpr sum_xy(env);
		for (IloInt j = 0; j < DEVICES; j++)
		{
			sum_xy.clear();
			for (IloInt i = 0; i < ACCESS; i++) {
				sum_xy += (x[i][j] + y[i][j]);
			}
			model.add(sum_xy >= 1);
		}

		//目标函数
		IloExpr obj(env);
		for (IloInt j = 0; j < DEVICES; j++) {
			for (IloInt i = 0; i < ACCESS; i++) {
				obj += Energy_cloud[i][j] * y[i][j];
			}
		}
		for (IloInt j = 0; j < DEVICES; j++) {
			for (IloInt i = 0; i < ACCESS; i++)
				obj += Energy_edge[i][j] * x[i][j];
		}
		model.add(IloMinimize(env, obj));

		//创建求解对象
		IloCplex cplex(model);
		if (!cplex.solve()) {
			cplex.exportModel("c:/cplex/cplexmodel.lp");
			env.error() << "Failed to optimize LP." << endl;
			throw(-1);
		}

		solution_value = cplex.getObjValue();

		printf("Solution value = %lf", solution_value);

		for (IloInt i = 0; i < ACCESS; i++)
		{
			IloNumArray xx(env);
			cplex.getValues(xx, x[i]);

			for (IloInt j = 0; j < DEVICES; j++) {
				xc[i][j] = xx[j];
			}

		}
		cout << endl;
		for (IloInt i = 0; i < ACCESS; i++)
		{
			IloNumArray yy(env);
			cplex.getValues(yy, y[i]);

			for (IloInt j = 0; j < DEVICES; j++) {
				yc[i][j] = yy[j];
			}
		}

	}
	catch (IloException& e) { cerr << "Concert exception caught:" << e << endl; }
	catch (...) { cerr << "Unknuwn exception caught" << endl; }
	env.end();

	// 更新结果
	cal_result(result, xc, yc, isPrint.p_xy);

}


//最小化云能耗
void cplex_cloud(Result* result)
{
	IloEnv env;
	double solution_value = 0;
	try {

		IloModel model(env);

		IloNumVarArray2 x(env, ACCESS);
		IloNumVarArray2 y(env, ACCESS);
		for (IloInt i = 0; i < ACCESS; i++)
		{
			x[i] = IloNumVarArray(env, DEVICES);
			y[i] = IloNumVarArray(env, DEVICES);
			for (IloInt j = 0; j < DEVICES; j++) {
				x[i][j] = IloNumVar(env, 0, 1, ILOINT);
				model.add(x[i][j]);

				y[i][j] = IloNumVar(env, 0, 1, ILOINT);
				model.add(y[i][j]);
			}
		}

		//约束8a
		IloExpr sum_edge_r(env);
		for (IloInt i = 0; i < ACCESS; i++)
		{
			for (IloInt k = 0; k < RESOURCE; k++)
			{
				sum_edge_r.clear();
				for (IloInt j = 0; j < DEVICES; j++)
					sum_edge_r += device[j].task.r[k] * x[i][j];
				model.add(sum_edge_r <= access[i].constrain_r[k]);
			}
		}


		//约束8b
		IloExpr sum_edge_BW(env);
		for (IloInt i = 0; i < ACCESS; i++)
		{
			sum_edge_BW.clear();
			for (IloInt j = 0; j < DEVICES; j++)
			{
				sum_edge_BW += device[j].task.r_BW * (x[i][j] + y[i][j]);
			}
			model.add(sum_edge_BW <= access[i].constrain_BW);
		}

		//约束8c
		IloExpr sum_cloud_r(env);
		for (IloInt k = 0; k < RESOURCE; k++)
		{
			sum_cloud_r.clear();
			for (IloInt i = 0; i < ACCESS; i++)
			{
				for (IloInt j = 0; j < DEVICES; j++)
				{
					sum_cloud_r += device[j].task.r[k] * y[i][j];
				}
			}
			model.add(sum_cloud_r <= cloud.constrain_r[k]);
		}

		//约束8d
		IloExpr sum_cloud_BW(env);
		for (IloInt i = 0; i < ACCESS; i++)
		{
			for (IloInt j = 0; j < DEVICES; j++) {
				sum_cloud_BW += device[j].task.r_BW * y[i][j];
			}
		}
		IloExpr fenzi(env);
		IloExpr fenmu(env);
		IloExpr temp(env);

		for (IloInt i = 0; i < ACCESS; i++)
		{
			temp.clear();
			for (IloInt j = 0; j < DEVICES; j++) {
				temp += device[j].task.r_BW * x[i][j];
			}
			try {
				sum_cloud_BW += access[i].sigma * temp;		
			}
			catch (exception e) {
				sum_cloud_BW += 0;
			}
		}
		model.add(sum_cloud_BW <= cloud.constrain_BW);

		//约束8e
		IloExpr sum_xy(env);
		for (IloInt j = 0; j < DEVICES; j++)
		{
			sum_xy.clear();
			for (IloInt i = 0; i < ACCESS; i++) {
				sum_xy += (x[i][j] + y[i][j]);
			}
			model.add(sum_xy >= 1);
		}

		//目标函数
		IloExpr obj(env);
		for (IloInt j = 0; j < DEVICES; j++) {
			for (IloInt i = 0; i < ACCESS; i++) {
				obj += Energy_cloud[i][j] * y[i][j];
			}
		}

		model.add(IloMinimize(env, obj));

		//创建求解对象
		IloCplex cplex(model);
		if (!cplex.solve()) {
			env.error() << "Failed to optimize LP." << endl;
			throw(-1);
		}

		solution_value = cplex.getObjValue();

		printf("Solution value = %lf", solution_value);

		cout << endl;
		for (IloInt i = 0; i < ACCESS; i++)
		{
			IloNumArray xx(env);
			cplex.getValues(xx, x[i]);

			for (IloInt j = 0; j < DEVICES; j++) {
				xc_cloud[i][j] = xx[j];
			}
		}
		cout << endl;
		for (IloInt i = 0; i < ACCESS; i++)
		{
			IloNumArray yy(env);
			cplex.getValues(yy, y[i]);
			for (IloInt j = 0; j < DEVICES; j++) {
				yc_cloud[i][j] = yy[j];
			}
		}

	}
	catch (IloException& e) { cerr << "Concert exception caught:" << e << endl; }
	catch (...) { cerr << "Unknuwn exception caught" << endl; }
	env.end();

	// 更新结果
	cal_result(result, xc_cloud, yc_cloud, isPrint.p_xy);

}



//最小化边能耗
void cplex_edge(Result* result)
{
	IloEnv env;
	double solution_value = 0;
	try {

		IloModel model(env);

		IloNumVarArray2 x(env, ACCESS);
		IloNumVarArray2 y(env, ACCESS);
		for (IloInt i = 0; i < ACCESS; i++)
		{
			x[i] = IloNumVarArray(env, DEVICES);
			y[i] = IloNumVarArray(env, DEVICES);
			for (IloInt j = 0; j < DEVICES; j++) {
				x[i][j] = IloNumVar(env, 0, 1, ILOINT);
				model.add(x[i][j]);

				y[i][j] = IloNumVar(env, 0, 1, ILOINT);
				model.add(y[i][j]);
			}
		}

		//约束8a
		IloExpr sum_edge_r(env);
		for (IloInt i = 0; i < ACCESS; i++)
		{
			for (IloInt k = 0; k < RESOURCE; k++)
			{
				sum_edge_r.clear();
				for (IloInt j = 0; j < DEVICES; j++)
					sum_edge_r += device[j].task.r[k] * x[i][j];
				model.add(sum_edge_r <= access[i].constrain_r[k]);
			}
		}


		//约束8b
		IloExpr sum_edge_BW(env);
		for (IloInt i = 0; i < ACCESS; i++)
		{
			sum_edge_BW.clear();
			for (IloInt j = 0; j < DEVICES; j++)
			{
				sum_edge_BW += device[j].task.r_BW * (x[i][j] + y[i][j]);
			}
			model.add(sum_edge_BW <= access[i].constrain_BW);
		}

		//约束8c
		IloExpr sum_cloud_r(env);
		for (IloInt k = 0; k < RESOURCE; k++)
		{
			sum_cloud_r.clear();
			for (IloInt i = 0; i < ACCESS; i++)
			{
				for (IloInt j = 0; j < DEVICES; j++)
				{
					sum_cloud_r += device[j].task.r[k] * y[i][j];
				}
			}
			model.add(sum_cloud_r <= cloud.constrain_r[k]);
		}

		//约束8d
		IloExpr sum_cloud_BW(env);
		for (IloInt i = 0; i < ACCESS; i++)
		{
			for (IloInt j = 0; j < DEVICES; j++) {
				sum_cloud_BW += device[j].task.r_BW * y[i][j];
			}
		}
		IloExpr fenzi(env);
		IloExpr fenmu(env);
		IloExpr temp(env);
		for (IloInt i = 0; i < ACCESS; i++)
		{
			temp.clear();
			for (IloInt j = 0; j < DEVICES; j++) {
				temp += device[j].task.r_BW * x[i][j];
			}
			try {
				sum_cloud_BW += access[i].sigma * temp;		
			}
			catch (exception e) {
				sum_cloud_BW += 0;
			}
		}
		model.add(sum_cloud_BW <= cloud.constrain_BW);

		//约束8e
		IloExpr sum_xy(env);
		for (IloInt j = 0; j < DEVICES; j++)
		{
			sum_xy.clear();
			for (IloInt i = 0; i < ACCESS; i++) {
				sum_xy += (x[i][j] + y[i][j]);
			}
			model.add(sum_xy >= 1);
		}

		//目标函数
		IloExpr obj(env);

		for (IloInt j = 0; j < DEVICES; j++) {
			for (IloInt i = 0; i < ACCESS; i++)
				obj += Energy_edge[i][j] * x[i][j];
		}
		model.add(IloMinimize(env, obj));

		//创建求解对象
		IloCplex cplex(model);
		if (!cplex.solve()) {
			env.error() << "Failed to optimize LP." << endl;
			throw(-1);
		}

		solution_value = cplex.getObjValue();

		printf("Solution value = %lf", solution_value);

		cout << endl;
		for (IloInt i = 0; i < ACCESS; i++)
		{
			IloNumArray xx(env);
			cplex.getValues(xx, x[i]);

			for (IloInt j = 0; j < DEVICES; j++) {
				xc_edge[i][j] = xx[j];
			}
		}
		cout << endl;
		for (IloInt i = 0; i < ACCESS; i++)
		{
			IloNumArray yy(env);
			cplex.getValues(yy, y[i]);

			for (IloInt j = 0; j < DEVICES; j++) {
				yc_edge[i][j] = yy[j];
			}
		}

	}
	catch (IloException& e) { cerr << "Concert exception caught:" << e << endl; }
	catch (...) { cerr << "Unknuwn exception caught" << endl; }
	env.end();

	// 更新结果
	cal_result(result, xc_edge, yc_edge, isPrint.p_xy);

}