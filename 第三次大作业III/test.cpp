#include<iostream>
#include<array>
#include"my_const.h"
#include"alpha_FPU.h"
#include"beta_FPU.h"
#include<fstream>
using namespace std;

//构造第一问的初始条件
//只有Q1非0,其他全为0
template<int n>
array<double, 2 * n> init_1(double Q1) {
	array<double, 2 * n> result;
	//从Q表象变换到q表象
	//坐标部分
	for (int i = 1; i <= n; i++) {
		result[i - 1] = Q1*sqrt(2.*(n+1)) / (n + 1)*sin(PI * 1 * i / double(n + 1));
	}
	for (int i = n; i < 2 * n; i++) {
		result[i] = 0;
	}
	return result;
}

//测试第二问
void test_2() {
	ofstream os;
	os.open("temp2.txt");
	//设置参数
	double alpha = 0.25;
	constexpr int n = 32;
	double delta_t = 0.6;//时间步长
	double Q1 = 4;
	//初始化
	alpha_FPU<n> FPU2(alpha, delta_t, init_1<n>(Q1));
	//一个周期T以及总时长tf
	double T = 2 * PI / FPU2.omega(1);
	double tf = 160 * 2 * PI / FPU2.omega(1);
	//迭代
	for (double t = 0; t < tf; t += delta_t) {
		//输出E1,E2,E3,E4
		os << t/T<<"\t"
			<< FPU2.Enegy(1) << "\t"
			<< FPU2.Enegy(2) << "\t"
			<< FPU2.Enegy(3) << "\t"
			<< FPU2.Enegy(4) << endl;
		//迭代一步
		FPU2.t_iter();
	}
}

//测试第4问
void test_4() {
	ofstream os;
	os.open("temp4.txt");
	//设置参数
	double alpha = 0.25;
	constexpr int n = 32;
	double delta_t = 0.1;//时间步长
	double Q1 = 4;
	//初始化
	alpha_FPU<n> FPU4(alpha, delta_t, init_1<n>(Q1));
	//一个周期T以及总时长tf
	double T = 2 * PI / FPU4.omega(1);
	double tf = 4000 * 2 * PI / FPU4.omega(1);
	//迭代
	for (double t = 0; t < tf; t += delta_t) {
		//输出E1,E2,E3,E4
		os << t / T << "\t"
			<< FPU4.Enegy(1) << endl;
		//迭代一步
		FPU4.t_iter();
	}
}

//第5问 E1的演化
void test_5_E1() {
	ofstream os;
	os.open("temp5_E1.txt");
	//设置参数
	double alpha = 0.25;
	constexpr int n = 32;
	double delta_t = 0.1;//时间步长
	double Q1 = 20;
	//初始化
	alpha_FPU<n> FPU4(alpha, delta_t, init_1<n>(Q1));
	//一个周期T以及总时长tf
	double T = 2 * PI / FPU4.omega(1);
	double tf = 4000 * 2 * PI / FPU4.omega(1);
	//迭代
	for (double t = 0; t < tf; t += delta_t) {
		//输出E1,E2,E3,E4
		os << t / T << "\t"
			<< FPU4.Enegy(1) << endl;
		//迭代一步
		FPU4.t_iter();
	}
}
//测试第5问平均能量E1,E2,E3,E4的演化
void test_5_average() {
	ofstream os;
	os.open("temp5_average.txt");
	//设置参数
	double alpha = 0.25;
	constexpr int n = 32;
	double delta_t = 0.1;//时间步长
	double Q1 = 20;
	//初始化
	alpha_FPU<n> FPU2(alpha, delta_t, init_1<n>(Q1));
	//一个周期T以及总时长tf
	double T = 2 * PI / FPU2.omega(1);
	double tf = 160 * 2 * PI / FPU2.omega(1);
	//迭代
	for (double t = 0; t < tf; t += delta_t) {
		//输出E1,E2,E3,E4
		os << t / T << "\t"
			<< FPU2.Enegy(1) << "\t"
			<< FPU2.Enegy(2) << "\t"
			<< FPU2.Enegy(3) << "\t"
			<< FPU2.Enegy(4) << endl;
		//迭代一步
		FPU2.t_iter();
	}
}

//第6问!!!!!!!!!!!!!!!!!!!!!!

//第7问
void test_7() {
	ofstream os;
	os.open("temp7.txt");
	//设置参数
	double beta = 0.25;
	constexpr int n = 32;
	double delta_t = 0.1;//时间步长
	double Q1 = 20;
	//初始化
	beta_FPU<n> FPU2(beta, delta_t, init_1<n>(Q1));
	//一个周期T以及总时长tf
	double T = 2 * PI / omega<n>(1);
	double tf = 160 * 2 * PI / omega<n>(1);
	//迭代
	for (double t = 0; t < tf; t += delta_t) {
		//输出E1,E2,E3,E4
		os << t / T << "\t"
			<< FPU2.Enegy(1) << "\t"
			<< FPU2.Enegy(2) << "\t"
			<< FPU2.Enegy(3) << "\t"
			<< FPU2.Enegy(4) << endl;
		//迭代一步
		FPU2.t_iter();
	}
}
int main() {
	test_2();
	system("pause");
}