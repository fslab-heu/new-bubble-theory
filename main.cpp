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
#include "global.h"             //��������
#include "bubble.h"            //��**********�����ݶ���************��


int main()
{
    #include "input.h"          //�����ļ�
    #include "protreat.h"      //��������
	while (t<Tend)
	{
		//ʱ�䲽
		allocator<double> alloc;
		double *Rmin = alloc.allocate(numb);
		for (int i = 0; i < numb; i++) { Rmin[i] =  bub[i].R; }
		dt = getMin(Rmin, numb)*2e-5;
		     
		induce_free(&bub[0],numb);    //������
		if (boundary == 1)  {induce_boun(&bub[0], numb, pos_boun, boun_norm);}  //�߽磨�����߽磩
		for (int i = 0; i < numb; i++)  {solve(&bub[i], dt);}  //������Ǩ�Ʒ���

       #include "output.h"               //����ļ�

		kk = kk + 1;
		t = t + dt;
		cout << kk << "\t" << "t = " << t << "\t" << endl;
	}

	return 0;
	
}

