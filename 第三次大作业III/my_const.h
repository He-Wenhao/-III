#pragma once

//用来统一储存一些常量
constexpr double E = 2.71828182845904523536028747;
constexpr double PI = 3.14159265358979323846264338;


//求k阶频率
template<int n>
double omega(int k) {
	return 2 * sin(PI*k / 2. / (n + 1));
}