FILE *fp3;
fp3 = fopen("thermodynamic.dat", "r");
fscanf(fp3, "%lf %s", &Rv, sstr);
fscanf(fp3, "%lf %s", &Rair, sstr);
//fscanf(fp3, "%lf %s", &kg, sstr);
fscanf(fp3, "%lf %s", &kq, sstr);
fscanf(fp3, "%lf %s", &T_ambient, sstr);
fscanf(fp3, "%lf %s", &M_air, sstr);
fscanf(fp3, "%lf %s", &M_vapour, sstr);
fscanf(fp3, "%lf %s", &NA, sstr);
fscanf(fp3, "%lf %s", &bounTi, sstr);
fscanf(fp3, "%lf %s", &bounTo, sstr);
fscanf(fp3, "%lf %s", &heatm, sstr);
fscanf(fp3, "%lf %s", &a_air, sstr);
fscanf(fp3, "%lf %s", &a_vapour, sstr);
fscanf(fp3, "%lf %s", &b_air, sstr);
fscanf(fp3, "%lf %s", &b_vapour, sstr);
fscanf(fp3, "%lf %s", &sigma_r, sstr);
fscanf(fp3, "%lf %s", &cv_air, sstr);
fscanf(fp3, "%lf %s", &cv_vapour, sstr);
fscanf(fp3, "%lf %s", &kv, sstr);
fscanf(fp3, "%lf %s", &kl, sstr);
fscanf(fp3, "%lf %s", &k_b, sstr);
fscanf(fp3, "%lf %s", &Dl, sstr);
fclose(fp3);

FILE *fp2;
fp2 = fopen("constants.dat", "r");
fscanf(fp2, "%lf %s", &gama, sstr);
fscanf(fp2, "%lf %s", &C, sstr);
fscanf(fp2, "%lf %s", &vis, sstr);
fscanf(fp2, "%lf %s", &sur, sstr);
fscanf(fp2, "%lf %s", &Pv, sstr);
fscanf(fp2, "%lf %s", &gg, sstr);
fscanf(fp2, "%lf %s", &rho, sstr);
fscanf(fp2, "%lf %s", &Cd, sstr);
fscanf(fp2, "%lf %s", &Ca, sstr);
fscanf(fp2, "%lf %s", &Pa, sstr);
fclose(fp2);
g[0] = 0; g[1] = 0; g[2] =-gg;

FILE *fp1;
fp1 = fopen("case.dat", "r");
fscanf(fp1, "%lf %s", &Tend, sstr);
fscanf(fp1, "%lf  %s", &sdt, sstr);
fscanf(fp1, "%lf,%lf,%lf  %s", &pos_measure[0], &pos_measure[1], &pos_measure[2], sstr);
fscanf(fp1, "%d  %s", &boundary, sstr);
fscanf(fp1, "%lf, %lf, %lf  %s", &pos_boun[0], &pos_boun[1], &pos_boun[2], sstr);
fscanf(fp1, "%lf, %lf, %lf  %s", &boun_norm[0], &boun_norm[1], &boun_norm[2], sstr);
fscanf(fp1, "%lf  %s", &alpha1, sstr);
fscanf(fp1, "%d  %s", &migration, sstr);
fscanf(fp1, "%d  %s", &numb, sstr);
//allocator<bubble> alloc;
//bubble *bub = alloc.allocate(numb);
vector<bubble> bub(numb);
for (int i = 0; i < numb; i++)
{
	fscanf(fp1, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf   %s", &arr[0], \
		&arr[1], &arr[2], &arr[3], &arr[4], &arr[5], &arr[6], &arr[7], &arr[8], &arr[9], sstr);
	bub[i] = bubble{ arr[0], arr[1], arr[2], arr[3], arr[4], arr[5], arr[6], arr[7], arr[8], arr[9] };
}
fclose(fp1); 



