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
		result[i - 1] = Q1*sqrt(2.*(n)) / (n + 1)*sin(PI * 1 * i / double(n + 1));
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
	double T = 2 * PI / omega<n>(1);
	double tf = 160 * 2 * PI / omega<n>(1);
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
	double T = 2 * PI / omega<n>(1);
	double tf = 4000 * 2 * PI / omega<n>(1);
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
	double T = 2 * PI / omega<n>(1);
	double tf = 4000 * 2 * PI / omega<n>(1);
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
	double delta_t = 0.01;//时间步长
	double Q1 = 20;
	//初始化
	alpha_FPU<n> FPU2(alpha, delta_t, init_1<n>(Q1));
	//平均值
	double ave1{ 0 }, ave2{ 0 }, ave3{ 0 }, ave4{ 0 };
	int n_iter = 0;//迭代次数
	//一个周期T以及总时长tf
	double T = 2 * PI / omega<n>(1);
	double tf = 4000 * 2 * PI / omega<n>(1);
	//迭代
	const int N_tol = tf / delta_t;
	const int delta_cout = N_tol / (2e4);
	for (double t = 0; t < tf; t += delta_t) {
		//更新平均值
		ave1 = (n_iter*ave1 + FPU2.Enegy(1)) / (n_iter + 1);
		ave2 = (n_iter*ave2 + FPU2.Enegy(2)) / (n_iter + 1);
		ave3 = (n_iter*ave3 + FPU2.Enegy(3)) / (n_iter + 1);
		ave4 = (n_iter*ave4 + FPU2.Enegy(4)) / (n_iter + 1);
		//输出E1,E2,E3,E4
		if (n_iter%delta_cout == 0) {
			os << t / T << "\t"
				<< ave1 << "\t"
				<< ave2 << "\t"
				<< ave3 << "\t"
				<< ave4 << endl;
			/*os << t / T << "\t"
				<< FPU2.Enegy(1) << "\t"
				<< FPU2.Enegy(2) << "\t"
				<< FPU2.Enegy(3) << "\t"
				<< FPU2.Enegy(4) << endl;*/
		}
		//迭代一步
		FPU2.t_iter();
		n_iter++;
	}
}

//
void testexpr() {
	ofstream os;
	os.open("temp.txt");
	//设置参数
	double alpha = 0.25;
	constexpr int n = 32;
	double delta_t = 0.01;//时间步长
	double Q1 = 20;
	//初始化
	alpha_FPU<n> FPU4(alpha, delta_t, init_1<n>(Q1));
	//一个周期T以及总时长tf
	double T = 2 * PI / omega<n>(1);
	double tf = 4000 * 2 * PI / omega<n>(1);
	//!!!!!!!!!!!
	const int N_tol = tf / delta_t;
	int n_iter = 0;
	const int delta_cout = N_tol / (2e4);
	double av1 = 0;
	//迭代
	for (double t = 0; t < tf; t += delta_t) {
		av1 = (n_iter*av1 + FPU4.Enegy(2)) / (n_iter + 1);
		//输出E1,E2,E3,E4
		if (n_iter%delta_cout == 0) {
			os << t / T << "\t"
				<< av1 << endl;
			//迭代一步
		}
			FPU4.t_iter();
		n_iter++;
	}
}


//测试第6问
void test_6_1() {
	ofstream os;
	os.open("temp6_1.txt");
	//设置参数
	double beta = 0.1;
	constexpr int n = 16;
	double delta_t = 0.1;//时间步长
	double Q1 = 40;
	//初始化
	beta_FPU<n> FPU6_1(beta, delta_t, init_1<n>(Q1));
	//一个周期T以及总时长tf
	double T = 2 * PI / omega<n>(1);
	double tf = 500 * 2 * PI / omega<n>(1);
	//迭代
	for (double t = 0; t < tf; t += delta_t) {
		//输出E1,E2,E3,E4
		os << t / T << "\t"
			<< FPU6_1.Enegy(1) << endl;
		//迭代一步
		FPU6_1.t_iter();
	}
}//测试第6问
void test_6_2() {
	ofstream os;
	os.open("temp6_2.txt");
	//设置参数
	double beta = 0.1;
	constexpr int n = 16;
	double delta_t = 0.1;//时间步长
	double Q1 = 5;
	//初始化
	beta_FPU<n> FPU6_1(beta, delta_t, init_1<n>(Q1));
	//一个周期T以及总时长tf
	double T = 2 * PI / omega<n>(1);
	double tf = 500 * 2 * PI / omega<n>(1);
	//迭代
	for (double t = 0; t < tf; t += delta_t) {
		//输出E1,E2,E3,E4
		os << t / T << "\t"
			<< FPU6_1.Enegy(1) << endl;
		//迭代一步
		FPU6_1.t_iter();
	}
}
int main() {
	//test_5_average();
	//testexpr();
	test_6_1();
	test_6_2();
	system("pause");
}