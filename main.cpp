#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>

using namespace std;
#include < string >
#include <direct.h>
#include<vector>
#include "math.h"
#include "global.h"             //参数定义
#include "bubble.h"            //（**********）气泡对象（************）


int main()
{
    #include "input.h"          //读入文件
    #include "protreat.h"      //储存数据
	while (t<Tend)
	{
		//时间步
		allocator<double> alloc;
		double *Rmin = alloc.allocate(numb);
		for (int i = 0; i < numb; i++) { Rmin[i] =  bub[i].R; }
		dt = getMin(Rmin, numb)*2e-5;
		     
		induce_free(&bub[0],numb);    //多气泡
		if (boundary == 1)  {induce_boun(&bub[0], numb, pos_boun, boun_norm);}  //边界（单个边界）
		for (int i = 0; i < numb; i++)  {solve(&bub[i], dt);}  //脉动、迁移方程

       #include "output.h"               //输出文件

		kk = kk + 1;
		t = t + dt;
		cout << kk << "\t" << "t = " << t << "\t" << endl;
	}

	return 0;
	
}

