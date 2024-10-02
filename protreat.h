vector<double> arrsa2;
for (int i = 0; i < numb; i++)
{
	bub[i].cal_adu();
	arrsa2.push_back(t);
	arrsa2.push_back(bub[i].R);
	arrsa2.push_back(bub[i].dR);
	arrsa2.push_back(bub[i].ddR);
	arrsa2.push_back(bub[i].HH);
	arrsa2.push_back(bub[i].dHH);
	arrsa2.push_back(bub[i].adu);
	bub[i].arrsa.push_back(arrsa2);
	arrsa2.clear();
}

char cwd[256];
string scwd;
string to_string(int val);
_getcwd(cwd, 256);
scwd = cwd;
vector<string>fpU(numb);
vector<ofstream>ofs(numb);
//ofstream ofs;
for (int i = 0; i < numb; i++)
{
	double nn = i + 1;
	fpU[i] = scwd + "\\" + "bubble" + to_string(nn) + ".dat";
	ofs[i].open(fpU[i], ios::out);
}

string fpP;
ofstream ofsP;
fpP = scwd + "\\" + "pressure" + ".dat";
ofsP.open(fpP, ios::out);
