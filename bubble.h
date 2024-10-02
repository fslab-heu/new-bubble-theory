class bubble
{
public:
	double ts0, T0, T, Percent, R, R0;
	double dR, ddR;  // 
	double u[3], du[3], ddu[3], adu;  // 
	double HH;   // 
	double dHH;   // 
	double Pref, rhog;   // 
	double P, dP, P0, P_ind, dP_ind, P_ind_out, dP_ind_out;   // 
	double gradp[3], pos[3];   // 
	double gradp_out[3];
	double alpha_m, percent0, percent_now, m, m0, Rg,Q,dm;
	double n_air, n_vapour, nt, T_liquid, TB;
	vector<vector<double>> arrsa;

	bubble();
	bubble(double t_start, double poss1, double poss2, double poss3, double R0s, double dR0s, double P0s, double T0s, double Percents, double alphas);
	double interpz(double arr[][2], int n, double xx);
	double cal_adu();
	double cal_Pref();
	double induceb(double fr[], int size);
	double thermoP();
	double thermoPb();
};
bubble::bubble() {} //初始化
bubble::bubble(double t_start, double poss1, double poss2, double poss3, double R0s, double dR0s, double P0s, double T0s, double Percents, double alphas)
{
	ts0 = t_start;
	pos[0] = poss1; pos[1] = poss2; pos[2] = poss3;
	Pref = Pa + rho * gg*abs(pos[2]);
	R0 = R0s;
	dR = dR0s;
	P0 = P0s;
	P = P0;
	ddR = 0; dP = 0;
	HH = (P - Pref) / rho;
	dHH = (0) / rho;
	R = R0; P_ind = 0; dP_ind = 0;
	P_ind_out = 0; dP_ind_out = 0;  adu = 0;
	for (int i = 0; i < 3; i++) { gradp[i] = 0; gradp_out[i] = 0; }
	for (int i = 0; i < 3; i++) { u[i] = 0; du[i] = 0; ddu[i] = 0; }
	percent0 = Percents; percent_now = percent0; alpha_m = alphas; T = T0s;
	Rg = Rv* percent0 + Rair*(1.0 - percent0);
	m0 = P0 * (1.3333333*3.1415926*R*R*R) / Rg / T; m = m0; dm = 0;
	rhog = m0/(1.3333333*3.1415926*R*R*R);
	n_air = m0*(1- percent0)/M_air*NA; 
	n_vapour = m0 *percent0 / M_vapour*NA;
	nt = n_air + n_vapour;
	TB = T;  T_liquid = 0.5*TB + 0.5*T_ambient;
	Q = kq*(T- T_ambient)*1.333333*3.1415926*R*R*R;

}
double bubble::induceb(double fr[], int size)   //气泡pos在流场fr中诱导的压力和速度
{
	P_ind_out = 0; dP_ind_out = 0;
	for (int i = 0; i < 3; i++)
	{
		gradp_out[i] = 0;
	}

	double dis[3];
	for (int i = 0; i < 3; i++)
	{
		dis[i] = pos[i] - fr[i];
	}
	double adis = sqrt(dis[0] * dis[0] + dis[1] * dis[1] + dis[2] * dis[2]);
	double tyc = t - abs(adis - R) / C; double Ryc = R; double dRyc = dR; double ddRyc = ddR; double HHyc = HH; double dHHyc = dHH; double aduyc = adu;
	if (tyc > 0 && kk > 5)
	{
		int ii = -1;
		for (int i = 0; i < kk-1; i++)
		{
			if (tyc >= arrsa[i][0] && tyc <= arrsa[i + 1][0]) {
				ii = i;
				break;
			}
		}
		if (ii < 0)
		{
			//Ryc = ((arrsa[kk - 2][0] - tyc)*arrsa[kk-3][1] - (arrsa[kk - 3][0] - tyc)*arrsa[kk - 2][1]) / (arrsa[kk-2][0] - arrsa[kk-3][0]);
			//dRyc = ((arrsa[kk - 2][0] - tyc)*arrsa[kk - 3][2] - (arrsa[kk - 3][0] - tyc)*arrsa[kk - 2][2]) / (arrsa[kk - 2][0] - arrsa[kk - 3][0]);
			//ddRyc = ((arrsa[kk - 2][0] - tyc)*arrsa[kk - 3][3] - (arrsa[kk - 3][0] - tyc)*arrsa[kk - 2][3]) / (arrsa[kk - 2][0] - arrsa[kk - 3][0]);
			//HHyc = ((arrsa[kk - 2][0] - tyc)*arrsa[kk - 3][4] - (arrsa[kk - 3][0] - tyc)*arrsa[kk - 2][4]) / (arrsa[kk - 2][0] - arrsa[kk - 3][0]);
			//dHHyc = ((arrsa[kk - 2][0] - tyc)*arrsa[kk - 3][5] - (arrsa[kk - 3][0] - tyc)*arrsa[kk - 2][5]) / (arrsa[kk - 2][0] - arrsa[kk - 3][0]);
			//aduyc = ((arrsa[kk - 2][0] - tyc)*arrsa[kk - 3][6] - (arrsa[kk - 3][0] - tyc)*arrsa[kk - 2][6]) / (arrsa[kk - 2][0] - arrsa[kk - 3][0]);
			/*cout << "Error in interp" << endl; exit(0);*/
		}
		else
		{
		Ryc = ((arrsa[ii + 1][0] - tyc)*arrsa[ii][1] - (arrsa[ii][0] - tyc)*arrsa[ii + 1][1]) / (arrsa[ii + 1][0] - arrsa[ii][0]);
		dRyc = ((arrsa[ii + 1][0] - tyc)*arrsa[ii][2] - (arrsa[ii][0] - tyc)*arrsa[ii + 1][2]) / (arrsa[ii + 1][0] - arrsa[ii][0]);
		ddRyc = ((arrsa[ii + 1][0] - tyc)*arrsa[ii][3] - (arrsa[ii][0] - tyc)*arrsa[ii + 1][3]) / (arrsa[ii + 1][0] - arrsa[ii][0]);
		HHyc = ((arrsa[ii + 1][0] - tyc)*arrsa[ii][4] - (arrsa[ii][0] - tyc)*arrsa[ii + 1][4]) / (arrsa[ii + 1][0] - arrsa[ii][0]);
		dHHyc = ((arrsa[ii + 1][0] - tyc)*arrsa[ii][5] - (arrsa[ii][0] - tyc)*arrsa[ii + 1][5]) / (arrsa[ii + 1][0] - arrsa[ii][0]);
		aduyc = ((arrsa[ii + 1][0] - tyc)*arrsa[ii][6] - (arrsa[ii][0] - tyc)*arrsa[ii + 1][6]) / (arrsa[ii + 1][0] - arrsa[ii][0]);
		}
	}
	P_ind_out = rho * Ryc / adis * (HHyc + 0.5*dRyc*dRyc+0.25*aduyc*aduyc);
	dP_ind_out = rho * (dRyc*HHyc + Ryc * dHHyc + 0.5*dRyc*dRyc*dRyc + Ryc*dRyc*ddRyc) / adis;
	for (int i = 0; i < 3; i++)
	{
		gradp_out[i] = -rho * Ryc / adis / adis * (HHyc + 0.5*dRyc*dRyc + 0.25*aduyc*aduyc)*dis[i] / adis;
	}
	return P_ind_out, dP_ind_out, gradp_out[3];
}
double bubble::thermoP()      //简化的热力学模型，待完善
{
	double Pvxx = exp(77.345 + 0.0057*T_liquid - 7235 / T_liquid) / pow(T_liquid, 8.2);
	double dm_eva = alpha_m * Pvxx / sqrt(2.0 * 3.1415926*Rv) / sqrt(T_liquid);
	double percent_now = (m0*percent0 + (m - m0)) / m;
	double dm_con = alpha_m * 1.0*P*percent_now / sqrt(2 * 3.1415926*Rv) / sqrt(abs(T));
	double dm = (dm_eva - dm_con)*(4.0 * 3.1415926*R * R);
	if (m + dm * dt<m0 && abs(m + dm * dt - m0) / m0 > percent0)    {dm = 0;}
	m = m + dm * dt;

	double cv = Rv / (gama - 1.0);
	double cp = cv + Rv;
	double dQ = kq * (T - T_liquid) * 4 * 3.1415926*R*R;
	//double Q = Q + dQ * dt;
	double dE = -4.0 * 3.1415926*R*R * dR*P - (dm*cv*T - dm * cp*T) - dQ;
	double dT = dE / m / cv;
	T = T + dT * dt;
	Rg = Rv * percent_now + Rair * (1 - percent_now);
	P = m * Rg*T / (4.0 / 3.0 * 3.1415926*R*R*R);
	return P; 
}
double bubble::thermoPb()    //论文中的热力学模型
{
	// 冷凝和蒸发速率
	double Pvxx = 133.0*exp(18.3036 - 3816.0/ (T_liquid-46.13)) ;
	double dm_eva = alpha_m * Pvxx / sqrt(2.0 * 3.1415926*Rv) / sqrt(T_liquid);
	double percent_now = n_vapour / nt;
	double dm_con = alpha_m * 1.0*P*percent_now / sqrt(2.0 * 3.1415926*Rv) / sqrt(abs(TB));
	 dm = (dm_eva - dm_con)*(4.0 * 3.1415926*R * R);
	if (m + dm * dt<m0 && abs(m + dm * dt - m0) / m0 > percent0) { dm = 0; }
	m = m + dm * dt;

	//气体分子数量
	double D = 1.76e-9; double KB = 6.737e9; double c0 = 4.8e23; double cs = 1000 * NA*(n_air*P / nt) / (M_vapour*KB);
	double dn_vapour = dm*NA / M_vapour;
	double dn_air = 4.0 * 3.1415926*R*D*(c0 - cs);
	n_vapour = n_vapour + dt * dn_vapour;
	n_air = n_air + dt * dn_air;
	nt = n_air + n_vapour;

	//气泡内部热力学边界层（可关）
	double a1 = 0.827; double alpha_e = 1.0; 
	//double kk = 0.02; double k_b = 1.38e-23; 
	double n1 = nt / (1.3333*3.1415926*R*R*R);
	double m_ba = (n_vapour*M_vapour + n_air * M_air) / NA / nt;
	double delta1 = 1.3333*3.1415926*R*R*R / sqrt(2.0) / 0.4e-18 / nt;
	double T1_r;
	if (TB != T)
	{
	    T1_r = (TB - T) / delta1 / bounTi;
		double kaf = -0.5 / k_b / n1 * sqrt(3.1415926*m_ba / 2.0 / k_b / TB) *(2 - a1 * alpha_e) / alpha_e * kv;
		TB = T_liquid + kaf * T1_r;
	}
	else
	{
		T1_r = 0;
	}
	

	//气泡外部热力学边界层（可关）
	//double kl = 0.5; double Dl = 1.43e-7; 
	double cv = Rv / (gama - 1); double cp = Rv / (1.25 - 1) + Rv; double Latent;
	if (T < 643.0)
	{
		Latent = 244281.0*pow((673.0 - 9.0 / 5.0 * (T - 273.0)) , 0.358);
	}	
	else
	{
		Latent = 244281.0*pow((673.0 - 9.0 / 5.0 * (643.0 - 273.0)) , 0.358);
	}
	double T1L_r = (kv*T1_r + dm/(4.0*3.1415926*R*R) * Latent) / kl;
	double e1L = 0.001*abs((T_liquid - TB) / T1L_r);
	double BL = T1L_r / 2.0 / e1L / (T_liquid - T_ambient);
	double CL = R + e1L;
	double AL = (T_liquid - T_ambient)*exp(BL*e1L* e1L);
	if (dR < 0)
	{
		if ((T_liquid - T_ambient)*T1L_r < 0)
		{
			double delta_e = -bounTo * (T_liquid - T_ambient) / T1L_r;
			double T1L_rex = (T_liquid - T_ambient)*exp(-delta_e / (T_ambient - T_liquid)*T1L_r)*(-T1L_r / (T_ambient - T_liquid));
			double xielv = (4.0*3.1415926*R*R * (-kl * T1L_r) - 4.0 * 3.1415926*pow((R + delta_e), 2.0) * (-kl * T1L_rex)) / (1.3333*3.1415926*rho*cp*(pow((R + delta_e), 3.0) - R * R*R));
			T_liquid = T_liquid + xielv * dt;
		}
		else if ((T_liquid - T_ambient)*T1L_r > 0)
		{
			double delta_e = bounTo * (e1L + 1 / sqrt(BL));
			double T1L_rex = AL * exp(-BL * pow((R + delta_e - CL), 2.0))*(-2.0 * BL*(R + delta_e - CL));
			double xielv = (4.0*3.1415926*R*R * (-kl * T1L_r) - 4.0 * 3.1415926*pow((R + delta_e), 2.0) * (-kl * T1L_rex)) / (1.3333*3.1415926*rho*cp*(pow((R + delta_e), 3.0) - R * R*R));
			T_liquid = T_liquid + xielv * dt;
		}
		else
		{
			T_liquid = T_liquid;
		}
	}
	else
   { T_liquid = 0.5*TB + 0.5*T_ambient; }
	if (bounTo==0)
	{
		T_liquid = T_ambient;
	}

	//气泡中心温度
	//double sigma_r = 5.67e-8; double cv_air = 1089.0; double cv_vapour = 717.0;
	double e_vapour_T = 6.272e-23*T - 5.247e-21;
	double e_vapour_Tl = 6.272e-23*T_liquid - 5.247e-21; double e_vapour_TB = 6.272e-23*TB - 5.247e-21; double e_air = 0;
	double dE = -4.0*3.1415926*R*R*dR*P + 4.0*3.1415926*R*R*(dm_eva*e_vapour_Tl - dm_con * e_vapour_TB)*NA / M_vapour \
		+ 4.0*3.1415926*R*R* 0.0*(T - T_liquid) + heatm *4.0*3.1415926*R*R* kv*T1_r + 0.0*sigma_r *4.0*3.1415926*R*R* (TB*TB*TB*TB - T*T*T*T);
	//double dE = -4.0*3.1415926*R*R*dR*P ;
	double cv1 = (n_air*cv_air + n_vapour * cv_vapour) / (n_air + n_vapour);
	 rhog = (n_air*M_air + n_vapour * M_vapour) / (1.3333*3.1415926*R*R*R) / NA;
	double dT = dE / (rhog *1.3333*3.1415926*R*R*R) / cv1;
	T = T + dT * dt;


	//气泡内压
	//double a_air = 0.1402; double a_vapour = 0.5536; double b_air = 3.753e-5; double b_vapour = 3.049e-5;
	double a = pow(sqrt(a_air)*n_air / nt+sqrt(a_vapour)*n_vapour / nt,2.0);
	double b = pow(sqrt(b_air)*n_air / nt+sqrt(b_vapour)*n_vapour / nt,2.0);
	double vv = NA * (1.3333*3.1415926*R*R*R) / nt;
	P = 8.314 * T / (vv - b) - a / vv /vv;

	//气泡内压时间导数
	double dnt = dn_air + dn_vapour;
	double da = 2.0 * (sqrt(a_air)*n_air / nt + sqrt(a_vapour)*n_vapour / nt)*(sqrt(a_air)*(dn_air*nt - n_air * dnt) / nt / nt + sqrt(a_vapour)*(dn_vapour*nt - n_vapour * dnt) / nt / nt);
	double db = 2.0 * (sqrt(b_air)*n_air / nt + sqrt(b_vapour)*n_vapour / nt)*(sqrt(b_air)*(dn_air*nt - n_air * dnt) / nt / nt + sqrt(b_vapour)*(dn_vapour*nt - n_vapour * dnt) / nt / nt);
	double dvv = NA * 4.0 * 3.1415926*R*R * dR / nt - NA * 1.3333*3.1415926*R*R*R * dnt / nt / nt;
	dP = 8.314 * dT / (vv - b) - 8.314 * T*(dvv - db) / (vv - b) / (vv - b) - da / vv / vv + a * 2.0 * vv*dvv / vv/vv/vv/vv;

	return P, dP, T, T_liquid,m,dm;
}
double bubble::cal_adu() //算矢量长度
{
	adu = sqrt(du[0]*du[0]+du[1]*du[1]+ du[2] * du[2]);
	return adu;
}
double bubble::cal_Pref() //气泡位置处静水压
{
	Pref = Pa + rho * gg*abs(pos[2]);
	return Pref;
}
double bubble::interpz(double arr[][2], int n, double xx)  //插值函数
{
	if (xx < arr[0][0])
	{
		cout << "Error in interp" << endl;
		cout << xx << endl;
		exit(0);
	}
	double resul;
	for (int i = 0; i < n - 1; i++)
	{
		if (xx >= arr[i][0] && xx <= arr[i + 1][0])
		{
			resul = ((arr[i + 1][0] - xx)*arr[i][1] - (arr[i][0] - xx)*arr[i + 1][1]) / (arr[i + 1][0] - arr[i][0]);
		}
		if (abs(xx - arr[i][0]) < 0.0001*abs(arr[i + 1][0] - arr[i][0]))
		{
			resul = arr[i][1];
		}
		if (abs(xx - arr[i + 1][0]) < 0.0001*abs(arr[i + 1][0] - arr[i][0]))
		{
			resul = arr[i + 1][1];
		}
	}
	return resul;
}


void induce_free(bubble *bub, int num)  //多气泡相互作用（仅同相，异相待完善）
{
	for (int i = 0; i < num; i++)
	{
		bub[i].P_ind = 0; bub[i].dP_ind = 0;
		for (int k = 0; k < 3; k++) { bub[i].gradp[k] = 0; }
		for (int j = 0; j < num; j++)
		{
			if (i != j)
			{
				double a[3] = { bub[i].pos[0], bub[i].pos[1], bub[i].pos[2] };
				bub[j].induceb(a, 3);
				bub[i].P_ind = bub[i].P_ind + bub[j].P_ind_out;
				bub[i].dP_ind = bub[i].dP_ind + bub[j].dP_ind_out;
				for (int k = 0; k < 3; k++) { bub[i].gradp[k] = bub[i].gradp[k] + bub[j].gradp_out[k]; }
			}
		}
	}
}

void induce_boun(bubble *bub, int num, double pos_boun[3], double boun_norm[3])   //边界效应（单个边界）
{
	for (int i = 0; i < num; i++)
	{
		for (int j = 0; j < num; j++)
		{
			double dot = (pos_boun[0] - bub[j].pos[0])*boun_norm[0] + (pos_boun[1] - bub[j].pos[1])*boun_norm[1] + (pos_boun[2] - bub[j].pos[2])*boun_norm[2];
			double pos_image[3] = { bub[j].pos[0] + 2.0*dot*boun_norm[0],  bub[j].pos[1] + 2.0*dot*boun_norm[1], bub[j].pos[2] + 2.0*dot*boun_norm[2] };
			bub[i].induceb(pos_image, 3);
			bub[i].P_ind = bub[i].P_ind + alpha1 * bub[i].P_ind_out;
			bub[i].dP_ind = bub[i].dP_ind + alpha1 * bub[i].dP_ind_out;
			for (int k = 0; k < 3; k++) { bub[i].gradp[k] = bub[i].gradp[k] - alpha1 * bub[i].gradp_out[k]; }
		}
	}
}

void solve(bubble *bub, double dt)   //脉动、迁移方程
{
	bub->cal_adu();
	bub->cal_Pref();
	double dotudu = bub->du[0]* bub->ddu[0] + bub->du[1] * bub->ddu[1] + bub->du[2] * bub->ddu[2];
	double dms = bub->dm / (4.0 * 3.1415926*bub->R*bub->R);
	double ddm = (bub->dm - bub->arrsa[kk][7]) / dt / (4.0 * 3.1415926*bub->R*bub->R);
	//bub->P =  bub->P0*pow((bub->R0 / bub->R), 3.0*gama) ;
	bub->thermoPb();
	bub->HH = 1.0 / rho * (bub->P - bub->P_ind - bub->Pref - 2.0*sur / bub->R - 4.0* vis / bub->R*(bub->dR - dms / rho)  );
	//bub->dHH = 1.0 / rho * (3.0*gama*bub->P0*pow((bub->R0 / bub->R), 3.0*gama - 1.0)*(-1.0)*bub->R0*pow(bub->R, -2.0)*bub->dR);
	//bub->dHH = (bub->HH - bub->arrsa[kk][4])/dt;
	bub->dHH = 1.0 / rho*(bub->dP - bub->dP_ind  + 2.0*sur / bub->R/ bub->R*bub->dR + 4.0*vis / bub->R / bub->R*bub->dR*(bub->dR - dms / rho) - 4.0* vis / bub->R*(bub->ddR - ddm/ rho)) ;

	bub->ddR = ((1.0 + bub->dR / C)*bub->HH + bub->R / C * bub->dHH + 0.25*(1.0+ bub->dR/C)*bub->adu*bub->adu + bub->R/2.0/C* dotudu - 1.5*(1.0 - bub->dR / 3.0 / C)*pow(bub->dR, 2.0)
		+ dms / rho*(bub->dR + dms / 2.0 / rho + bub->dR * dms / 2.0 / rho / C) + ddm* bub->R / rho*(1 - bub->dR / C + dms / rho / C) )
		/ (1.0 - bub->dR / C + dms / rho / C) / bub->R;

	bub->dR = bub->dR + dt * bub->ddR;
	bub->R = bub->R + dt * bub->dR;

	if (migration == 1)
	{
		for (int i = 0; i < 3; i++)
		{
			double adu1 = sqrt(bub->du[0] * bub->du[0] + bub->du[1] * bub->du[1] + bub->du[2] * bub->du[2]);
			bub->ddu[i] = -bub->R / rho * (-bub->gradp[i] + g[i] * rho) - 3.0 / 8.0*Cd*bub->du[i] * adu1 - (3.0*Ca*bub->dR + 0 * bub->R)*bub->du[i] 
				+bub->rhog / rho * bub->R*g[i] - 3.0*dms/rho* bub->du[i];
			bub->ddu[i] = bub->ddu[i] / (Ca*bub->R + bub->rhog /rho* bub->R);
			bub->du[i] = bub->du[i] + dt * bub->ddu[i];
			bub->u[i] = bub->u[i] + dt * bub->du[i];
			bub->pos[i] = bub->pos[i] + dt * bub->du[i];
		}
	}
}