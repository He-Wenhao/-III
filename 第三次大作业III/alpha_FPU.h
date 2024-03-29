﻿#pragma once
#include<array>
#include"my_const.h"
using namespace std;
//n代表自由度个数
template<int n>
struct alpha_FPU {
	double alpha;//三次非谐系数
	double delta_t;//迭代步长
	array<double, 2 * n> xp;//相空间坐标(x,p)
	//初始化
	alpha_FPU<n>(double _alpha, double _delta_t, array<double ,2 * n> _xp) {
		alpha = _alpha; 
		delta_t = _delta_t;
		xp = _xp;
	}
	array<double, 2 * n> diff_xp(array<double, 2 * n> xp);
	void t_iter();//完成一次迭代
	double Enegy(int k);//求第k阶能量
};




//重载array的+,*运算符
template<typename T,int K>
array<T, K> operator+(array<T, K> a, array<T, K> b) {
	array<T, K> result;
	for (int i = 0; i < K; i++) {
		result.at(i) = a.at(i) + b.at(i);
	}
	return result;
}

template<typename T, int K>
array<T, K> operator*(T a, array<T, K> b) {
	array<T, K> result;
	for (int i = 0; i < K; i++) {
		result.at(i) = a*b.at(i);
	}
	return result;
}





//利用哈密顿方程求dq/dt,dp/dt
template<int n>
array<double, 2 * n> alpha_FPU<n>::diff_xp(array<double, 2 * n> xp0) {
	array<double, 2 * n> result;
	//dq/dt=dH/dp
	for (int i = 0; i < n; i++) {
		result[i] = xp0[n + i];
	}
	double deltai, deltai1;
	//dp/dt=-dH/dq
	for (int i = n; i < 2 * n; i++) {
		if (i == n) {
			deltai= xp0.at(i - n) - xp0.at(i + 1 - n);
			deltai1 = xp0.at(i - n);
		}
		else if(i==2*n-1){
			deltai = xp0.at(n - 1);
			deltai1 = xp0.at(i - n) - xp0.at(i - 1 - n);
		}
		else {
			deltai = xp0[i - n] - xp0[i + 1 - n];
			deltai1 = xp0[i - n] - xp0.at(i - 1 - n);
		}
		result[i] = -deltai - deltai1 - alpha * pow(deltai, 2) + alpha * pow(deltai1, 2);
	}
	return result;
}

//利用四阶龙格库塔法做一次迭代
template<int n>
void alpha_FPU<n>::t_iter() {
	using a2n = array<double, 2 * n>;
	a2n k1 = diff_xp(xp);
	a2n k2 = diff_xp(xp + delta_t / 2.*k1);
	a2n k3 = diff_xp(xp + delta_t / 2.*k2);
	a2n k4 = diff_xp(xp + delta_t * k3);
	xp = xp + delta_t / 6.*(k1 + 2.*k2 + 2.*k3 + k4);
}


//求k阶能量
template<int n>
double alpha_FPU<n>::Enegy(int k) {
	//广义坐标坐标Qk
	double Qk = 0.;
	for (int j = 1; j < n+1; j++) {
		Qk += sin(PI*k*j / double(n+1 ))*xp[j - 1];
	}
	Qk = sqrt(2. / double(n))*Qk;
	//广义速度Qk_
	double Qk_ = 0;
	for (int j = 1; j < n + 1; j++) {
		Qk_ += sin(PI*k*j / double(n+1))*xp[j - 1 + n];
	}
	Qk_ = sqrt(2. / double(n))*Qk_;
	return 0.5*pow(Qk_, 2) + 0.5*pow(Qk*omega<n>(k), 2);
}