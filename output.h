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
	arrsa2.push_back(bub[i].dm);
	bub[i].arrsa.push_back(arrsa2);
	arrsa2.clear();

	ofs[i] << t << "," << bub[i].R << "," << bub[i].u[2] << endl;
}

P_loc = 0;
for (int i = 0; i < numb; i++)
{
	bub[i].cal_adu();
	bub[i].induceb(pos_measure, 3);
	P_loc = P_loc + bub[i].P_ind_out;
}

if(boundary == 1)
{
     for (int i = 0; i < numb; i++)
      {
             double dot = (pos_boun[0] - bub[i].pos[0])*boun_norm[0] + (pos_boun[1] - bub[i].pos[1])*boun_norm[1] + (pos_boun[2] - bub[i].pos[2])*boun_norm[2];
             double pos_image[3] = { bub[i].pos[0] + 2.0*dot*boun_norm[0],  bub[i].pos[1] + 2.0*dot*boun_norm[1], bub[i].pos[2] + 2.0*dot*boun_norm[2] };
			 for (int jj = 0; jj < 3; jj++) { pos_measureSave[jj] = bub[i].pos[jj]; }
			 for (int jj = 0; jj < 3; jj++) { bub[i].pos[jj] = pos_image[jj]; }
			 bub[i].induceb(pos_measure, 3);
			 P_loc = P_loc + alpha1* bub[i].P_ind_out;
			 for (int jj = 0; jj < 3; jj++) { bub[i].pos[jj] = pos_measureSave[jj]; }
     }
}

ofsP << t << "," << P_loc << endl;
