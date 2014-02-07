# SUSY Les Houches Accord 2.beta - MSSM spectrum + Decays
# SPheno module generated by SARAH
# ----------------------------------------------------------------------
# SPheno v3.1.11 
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# W. Porod, F.Staub, Comput.Phys.Commun.183 (2012) 2458-2469, arXiv:1104.1573
# SARAH: SARAHVER
# F. Staub, arXiv:0806.0538 (online manual)
# F. Staub, Comput. Phys. Commun. 181 (2010) 1077-1086, arXiv:0909.2863
# F. Staub, Comput. Phys. Commun. 182 (2011)  808-833, arXiv:1002.0840
# F. Staub, arXiv:1207.0906 
#  
# in case of problems send email to porod@physik.uni-wuerzburg.de 
# or to fnstaub@th.physik.uni-bonn.de 
# ----------------------------------------------------------------------
# Created: 23.05.2013,  17:04
Block SPINFO         # Program information
     1   SPhenoSARAH      # spectrum calculator
     2   v3.1.11     # version number of SPheno
     4   SARAHVER    # version number of SARAH
Block MODSEL  # Input parameters
     1  1   #  GUT scale input
     2  1   #  Boundary conditions 
Block MINPAR  # Input parameters
    1    2.40000000E+02  # m0
    2    9.70000000E+02  # m12
    3    1.60000000E+01  # TanBeta
    4    1.00000000E+00  # SignumMu
    5   -1.86000000E+03  # Azero
Block gaugeGUT Q=  1.34632464E+16  # (GUT scale)
   1    7.04308381E-01  # g1(Q)^DRbar
   2    7.04308381E-01  # g2(Q)^DRbar
   3    7.01248136E-01  # g3(Q)^DRbar
Block SMINPUTS  # SM parameters
         1     1.27931489E+02  # alpha_em^-1(MZ)^MSbar
         2     1.16639000E-05  # G_mu [GeV^-2]
         3     1.19000000E-01  # alpha_s(MZ)^MSbar
         4     9.11876000E+01  # m_Z(pole)
         5     4.20000000E+00  # m_b(m_b), MSbar
         6     1.73200000E+02  # m_t(pole)
         7     1.77700000E+00  # m_tau(pole)
Block GAUGE Q=  1.00000000E+03  # (SUSY scale)
   1    3.61967819E-01  # g1
   2    6.37320111E-01  # g2
   3    1.04480072E+00  # g3
Block HMIX Q=  1.00000000E+03  # (SUSY scale)
   1    2.59213754E+03  # Mu
 101    1.22878040E+05  # Bmu
 102    1.52951498E+01  # vd
 103    2.44722397E+02  # vu
   3    2.45199904E+02  # v
  10    1.50836574E+00  # betaH
  11    3.07772927E+00  # alphaH
Block MSOFT Q=  1.00000000E+03  # (SUSY scale)
  21   -6.20991922E+06  # mHd2
  22   -6.77643905E+06  # mHu2
   1    4.17978616E+02  # M1
   2    7.65447368E+02  # M2
   3    2.11437357E+03  # M3
Block PHASES Q=  1.00000000E+03  # (SUSY scale)
   1    1.00000000E+00  # pG
Block TREEHMIX Q=  1.00000000E+03  # (SUSY scale)
   1    2.60278756E+03  # Mu
 101    3.50286657E+04  # Bmu
Block Yu Q=  1.00000000E+03  # (SUSY Scale)
  1  1     8.31985232E-06   # Real(Yu(1,1),dp)
  1  2     0.00000000E+00   # Real(Yu(1,2),dp)
  1  3     0.00000000E+00   # Real(Yu(1,3),dp)
  2  1     0.00000000E+00   # Real(Yu(2,1),dp)
  2  2     3.52207090E-03   # Real(Yu(2,2),dp)
  2  3     0.00000000E+00   # Real(Yu(2,3),dp)
  3  1     0.00000000E+00   # Real(Yu(3,1),dp)
  3  2     0.00000000E+00   # Real(Yu(3,2),dp)
  3  3     8.48865534E-01   # Real(Yu(3,3),dp)
Block Yd Q=  1.00000000E+03  # (SUSY Scale)
  1  1     1.87906624E-04   # Real(Yd(1,1),dp)
  1  2     0.00000000E+00   # Real(Yd(1,2),dp)
  1  3     0.00000000E+00   # Real(Yd(1,3),dp)
  2  1     0.00000000E+00   # Real(Yd(2,1),dp)
  2  2     3.94603731E-03   # Real(Yd(2,2),dp)
  2  3     0.00000000E+00   # Real(Yd(2,3),dp)
  3  1     0.00000000E+00   # Real(Yd(3,1),dp)
  3  2     0.00000000E+00   # Real(Yd(3,2),dp)
  3  3     1.93909781E-01   # Real(Yd(3,3),dp)
Block Ye Q=  1.00000000E+03  # (SUSY Scale)
  1  1     4.73258869E-05   # Real(Ye(1,1),dp)
  1  2     0.00000000E+00   # Real(Ye(1,2),dp)
  1  3     0.00000000E+00   # Real(Ye(1,3),dp)
  2  1     0.00000000E+00   # Real(Ye(2,1),dp)
  2  2     9.78546846E-03   # Real(Ye(2,2),dp)
  2  3     0.00000000E+00   # Real(Ye(2,3),dp)
  3  1     0.00000000E+00   # Real(Ye(3,1),dp)
  3  2     0.00000000E+00   # Real(Ye(3,2),dp)
  3  3     1.64638768E-01   # Real(Ye(3,3),dp)
Block Tu Q=  1.00000000E+03  # (SUSY Scale)
  1  1    -2.81877807E-02   # Real(Tu(1,1),dp)
  1  2     0.00000000E+00   # Real(Tu(1,2),dp)
  1  3     0.00000000E+00   # Real(Tu(1,3),dp)
  2  1     0.00000000E+00   # Real(Tu(2,1),dp)
  2  2    -1.19327349E+01   # Real(Tu(2,2),dp)
  2  3     0.00000000E+00   # Real(Tu(2,3),dp)
  3  1     0.00000000E+00   # Real(Tu(3,1),dp)
  3  2     0.00000000E+00   # Real(Tu(3,2),dp)
  3  3    -1.95716308E+03   # Real(Tu(3,3),dp)
Block Td Q=  1.00000000E+03  # (SUSY Scale)
  1  1    -8.26202604E-01   # Real(Td(1,1),dp)
  1  2     0.00000000E+00   # Real(Td(1,2),dp)
  1  3     0.00000000E+00   # Real(Td(1,3),dp)
  2  1     0.00000000E+00   # Real(Td(2,1),dp)
  2  2    -1.73501295E+01   # Real(Td(2,2),dp)
  2  3     0.00000000E+00   # Real(Td(2,3),dp)
  3  1     0.00000000E+00   # Real(Td(3,1),dp)
  3  2     0.00000000E+00   # Real(Td(3,2),dp)
  3  3    -7.73536072E+02   # Real(Td(3,3),dp)
Block Te Q=  1.00000000E+03  # (SUSY Scale)
  1  1    -1.11873451E-01   # Real(Te(1,1),dp)
  1  2     0.00000000E+00   # Real(Te(1,2),dp)
  1  3     0.00000000E+00   # Real(Te(1,3),dp)
  2  1     0.00000000E+00   # Real(Te(2,1),dp)
  2  2    -2.31300768E+01   # Real(Te(2,2),dp)
  2  3     0.00000000E+00   # Real(Te(2,3),dp)
  3  1     0.00000000E+00   # Real(Te(3,1),dp)
  3  2     0.00000000E+00   # Real(Te(3,2),dp)
  3  3    -3.80865825E+02   # Real(Te(3,3),dp)
Block MSQ2 Q=  1.00000000E+03  # (SUSY Scale)
  1  1     3.70394342E+06   # Real(mq2(1,1),dp)
  1  2     0.00000000E+00   # Real(mq2(1,2),dp)
  1  3     0.00000000E+00   # Real(mq2(1,3),dp)
  2  1     0.00000000E+00   # Real(mq2(2,1),dp)
  2  2     3.70391355E+06   # Real(mq2(2,2),dp)
  2  3     0.00000000E+00   # Real(mq2(2,3),dp)
  3  1     0.00000000E+00   # Real(mq2(3,1),dp)
  3  2     0.00000000E+00   # Real(mq2(3,2),dp)
  3  3     3.50437649E+06   # Real(mq2(3,3),dp)
Block MSL2 Q=  1.00000000E+03  # (SUSY Scale)
  1  1     4.74286928E+05   # Real(ml2(1,1),dp)
  1  2     0.00000000E+00   # Real(ml2(1,2),dp)
  1  3     0.00000000E+00   # Real(ml2(1,3),dp)
  2  1     0.00000000E+00   # Real(ml2(2,1),dp)
  2  2     4.74313634E+05   # Real(ml2(2,2),dp)
  2  3     0.00000000E+00   # Real(ml2(2,3),dp)
  3  1     0.00000000E+00   # Real(ml2(3,1),dp)
  3  2     0.00000000E+00   # Real(ml2(3,2),dp)
  3  3     4.82709460E+05   # Real(ml2(3,3),dp)
Block MSD2 Q=  1.00000000E+03  # (SUSY Scale)
  1  1     3.37969275E+06   # Real(md2(1,1),dp)
  1  2     0.00000000E+00   # Real(md2(1,2),dp)
  1  3     0.00000000E+00   # Real(md2(1,3),dp)
  2  1     0.00000000E+00   # Real(md2(2,1),dp)
  2  2     3.37965159E+06   # Real(md2(2,2),dp)
  2  3     0.00000000E+00   # Real(md2(2,3),dp)
  3  1     0.00000000E+00   # Real(md2(3,1),dp)
  3  2     0.00000000E+00   # Real(md2(3,2),dp)
  3  3     3.29655843E+06   # Real(md2(3,3),dp)
Block MSU2 Q=  1.00000000E+03  # (SUSY Scale)
  1  1     3.41904961E+06   # Real(mu2(1,1),dp)
  1  2     0.00000000E+00   # Real(mu2(1,2),dp)
  1  3     0.00000000E+00   # Real(mu2(1,3),dp)
  2  1     0.00000000E+00   # Real(mu2(2,1),dp)
  2  2     3.41903013E+06   # Real(mu2(2,2),dp)
  2  3     0.00000000E+00   # Real(mu2(2,3),dp)
  3  1     0.00000000E+00   # Real(mu2(3,1),dp)
  3  2     0.00000000E+00   # Real(mu2(3,2),dp)
  3  3     3.07567278E+06   # Real(mu2(3,3),dp)
Block MSE2 Q=  1.00000000E+03  # (SUSY Scale)
  1  1     1.81204965E+05   # Real(me2(1,1),dp)
  1  2     0.00000000E+00   # Real(me2(1,2),dp)
  1  3     0.00000000E+00   # Real(me2(1,3),dp)
  2  1     0.00000000E+00   # Real(me2(2,1),dp)
  2  2     1.81256920E+05   # Real(me2(2,2),dp)
  2  3     0.00000000E+00   # Real(me2(2,3),dp)
  3  1     0.00000000E+00   # Real(me2(3,1),dp)
  3  2     0.00000000E+00   # Real(me2(3,2),dp)
  3  3     1.97628473E+05   # Real(me2(3,3),dp)
Block MASS  # Mass spectrum
#   PDG code      mass          particle
   1000001     1.82584627E+03  # Sd_1
   1000003     1.86168704E+03  # Sd_2
   1000005     1.86170171E+03  # Sd_3
   2000001     1.89011950E+03  # Sd_4
   2000003     1.95244705E+03  # Sd_5
   2000005     1.95245287E+03  # Sd_6
   1000002     1.69829043E+03  # Su_1
   1000004     1.87175465E+03  # Su_2
   1000006     1.87176468E+03  # Su_3
   2000002     1.93454365E+03  # Su_4
   2000004     1.95101707E+03  # Su_5
   2000006     1.95102147E+03  # Su_6
   1000011     4.17565372E+02  # Se_1
   1000013     4.33385575E+02  # Se_2
   1000015     4.33452095E+02  # Se_3
   2000011     6.96993378E+02  # Se_4
   2000013     6.97043465E+02  # Se_5
   2000015     7.11049852E+02  # Se_6
   1000012     6.92240767E+02  # Sv_1
   1000014     6.92241866E+02  # Sv_2
   1000016     6.92703833E+02  # Sv_3
        25     1.18244153E+02  # hh_1
        35     7.60541958E+02  # hh_2
        36     7.60784365E+02  # Ah_2
        37     7.65180509E+02  # Hpm_2
        23     9.11876000E+01  # VZ
        24     8.03850000E+01  # VWm
         1     5.00000000E-03  # Fd_1
         3     1.05000000E-01  # Fd_2
         5     4.20000000E+00  # Fd_3
         2     3.00000000E-03  # Fu_1
         4     1.27000000E+00  # Fu_2
         6     1.73200000E+02  # Fu_3
        11     5.10998910E-04  # Fe_1
        13     1.05658000E-01  # Fe_2
        15     1.77700000E+00  # Fe_3
   1000021     2.14061152E+03  # Glu
   1000022     4.15776159E+02  # Chi_1
   1000023     7.90390510E+02  # Chi_2
   1000025     2.57889555E+03  # Chi_3
   1000035     2.58037030E+03  # Chi_4
   1000024     7.90591661E+02  # Cha_1
   1000037     2.58083147E+03  # Cha_2
Block LSP  # LSP and NLSP
  1   1000022   # LSP 
  2   1000011   # NLSP 
Block FITTINOCHECK  # 
   1    1.00000000E+00  # GUT-scale mHd2 + mu^2 > 0 
   2    1.00000000E+00  # GUT-scale mHu2 + mu^2 > 0 
   3    0.00000000E+00  # GUT-scale ( mHd2 + mu^2 )(mHu2 + mu^2) - Bmu^2 > 0 
Block DSQMIX Q=  1.00000000E+03  # ()
  1  1    -0.00000000E+00   # Real(ZD(1,1),dp)
  1  2    -0.00000000E+00   # Real(ZD(1,2),dp)
  1  3    -4.20925619E-01   # Real(ZD(1,3),dp)
  1  4    -0.00000000E+00   # Real(ZD(1,4),dp)
  1  5    -0.00000000E+00   # Real(ZD(1,5),dp)
  1  6    -9.07095157E-01   # Real(ZD(1,6),dp)
  2  1    -0.00000000E+00   # Real(ZD(2,1),dp)
  2  2    -5.73606213E-03   # Real(ZD(2,2),dp)
  2  3    -0.00000000E+00   # Real(ZD(2,3),dp)
  2  4    -0.00000000E+00   # Real(ZD(2,4),dp)
  2  5    -9.99983549E-01   # Real(ZD(2,5),dp)
  2  6    -0.00000000E+00   # Real(ZD(2,6),dp)
  3  1    -2.73168702E-04   # Real(ZD(3,1),dp)
  3  2    -0.00000000E+00   # Real(ZD(3,2),dp)
  3  3    -0.00000000E+00   # Real(ZD(3,3),dp)
  3  4    -9.99999963E-01   # Real(ZD(3,4),dp)
  3  5    -0.00000000E+00   # Real(ZD(3,5),dp)
  3  6    -0.00000000E+00   # Real(ZD(3,6),dp)
  4  1     0.00000000E+00   # Real(ZD(4,1),dp)
  4  2     0.00000000E+00   # Real(ZD(4,2),dp)
  4  3    -9.07095157E-01   # Real(ZD(4,3),dp)
  4  4     0.00000000E+00   # Real(ZD(4,4),dp)
  4  5     0.00000000E+00   # Real(ZD(4,5),dp)
  4  6     4.20925619E-01   # Real(ZD(4,6),dp)
  5  1     0.00000000E+00   # Real(ZD(5,1),dp)
  5  2    -9.99983549E-01   # Real(ZD(5,2),dp)
  5  3     0.00000000E+00   # Real(ZD(5,3),dp)
  5  4     0.00000000E+00   # Real(ZD(5,4),dp)
  5  5     5.73606213E-03   # Real(ZD(5,5),dp)
  5  6     0.00000000E+00   # Real(ZD(5,6),dp)
  6  1    -9.99999963E-01   # Real(ZD(6,1),dp)
  6  2     0.00000000E+00   # Real(ZD(6,2),dp)
  6  3     0.00000000E+00   # Real(ZD(6,3),dp)
  6  4     2.73168702E-04   # Real(ZD(6,4),dp)
  6  5     0.00000000E+00   # Real(ZD(6,5),dp)
  6  6     0.00000000E+00   # Real(ZD(6,6),dp)
Block SNUMIX Q=  1.00000000E+03  # ()
  1  1     0.00000000E+00   # Real(ZV(1,1),dp)
  1  2     1.00000000E+00   # Real(ZV(1,2),dp)
  1  3     0.00000000E+00   # Real(ZV(1,3),dp)
  2  1     1.00000000E+00   # Real(ZV(2,1),dp)
  2  2     0.00000000E+00   # Real(ZV(2,2),dp)
  2  3     0.00000000E+00   # Real(ZV(2,3),dp)
  3  1     0.00000000E+00   # Real(ZV(3,1),dp)
  3  2     0.00000000E+00   # Real(ZV(3,2),dp)
  3  3     1.00000000E+00   # Real(ZV(3,3),dp)
Block USQMIX Q=  1.00000000E+03  # ()
  1  1     0.00000000E+00   # Real(ZU(1,1),dp)
  1  2     9.99110178E-19   # Real(ZU(1,2),dp)
  1  3    -4.80464320E-01   # Real(ZU(1,3),dp)
  1  4     2.27171353E-31   # Real(ZU(1,4),dp)
  1  5     4.04275131E-16   # Real(ZU(1,5),dp)
  1  6    -8.77014274E-01   # Real(ZU(1,6),dp)
  2  1     2.37195035E-16   # Real(ZU(2,1),dp)
  2  2     7.44403965E-03   # Real(ZU(2,2),dp)
  2  3    -1.30503867E-15   # Real(ZU(2,3),dp)
  2  4     1.34899598E-11   # Real(ZU(2,4),dp)
  2  5     9.99972293E-01   # Real(ZU(2,5),dp)
  2  6     1.17591687E-15   # Real(ZU(2,6),dp)
  3  1     1.75852998E-05   # Real(ZU(3,1),dp)
  3  2    -9.90058681E-14   # Real(ZU(3,2),dp)
  3  3     1.75855107E-26   # Real(ZU(3,3),dp)
  3  4     1.00000000E+00   # Real(ZU(3,4),dp)
  3  5    -1.34895966E-11   # Real(ZU(3,5),dp)
  3  6    -1.58521828E-26   # Real(ZU(3,6),dp)
  4  1    -0.00000000E+00   # Real(ZU(4,1),dp)
  4  2    -9.73438139E-17   # Real(ZU(4,2),dp)
  4  3    -8.77014274E-01   # Real(ZU(4,3),dp)
  4  4    -2.21334207E-29   # Real(ZU(4,4),dp)
  4  5    -1.70884636E-15   # Real(ZU(4,5),dp)
  4  6     4.80464320E-01   # Real(ZU(4,6),dp)
  5  1    -4.78970738E-16   # Real(ZU(5,1),dp)
  5  2     9.99972293E-01   # Real(ZU(5,2),dp)
  5  3    -7.51792009E-17   # Real(ZU(5,3),dp)
  5  4    -1.41395837E-15   # Real(ZU(5,4),dp)
  5  5    -7.44403965E-03   # Real(ZU(5,5),dp)
  5  6     3.88939690E-17   # Real(ZU(5,6),dp)
  6  1    -1.00000000E+00   # Real(ZU(6,1),dp)
  6  2    -4.78933769E-16   # Real(ZU(6,2),dp)
  6  3     3.57064226E-32   # Real(ZU(6,3),dp)
  6  4     1.75852998E-05   # Real(ZU(6,4),dp)
  6  5     3.53534043E-18   # Real(ZU(6,5),dp)
  6  6    -1.84728163E-32   # Real(ZU(6,6),dp)
Block SELMIX Q=  1.00000000E+03  # ()
  1  1    -0.00000000E+00   # Real(ZE(1,1),dp)
  1  2    -0.00000000E+00   # Real(ZE(1,2),dp)
  1  3    -2.41157512E-01   # Real(ZE(1,3),dp)
  1  4    -0.00000000E+00   # Real(ZE(1,4),dp)
  1  5    -0.00000000E+00   # Real(ZE(1,5),dp)
  1  6    -9.70485989E-01   # Real(ZE(1,6),dp)
  2  1    -1.97219094E-16   # Real(ZE(2,1),dp)
  2  2    -1.55186263E-02   # Real(ZE(2,2),dp)
  2  3    -0.00000000E+00   # Real(ZE(2,3),dp)
  2  4    -7.64186618E-17   # Real(ZE(2,4),dp)
  2  5    -9.99879579E-01   # Real(ZE(2,5),dp)
  2  6    -0.00000000E+00   # Real(ZE(2,6),dp)
  3  1     7.50775896E-05   # Real(ZE(3,1),dp)
  3  2    -1.18637219E-18   # Real(ZE(3,2),dp)
  3  3     0.00000000E+00   # Real(ZE(3,3),dp)
  3  4     9.99999997E-01   # Real(ZE(3,4),dp)
  3  5    -7.64242605E-17   # Real(ZE(3,5),dp)
  3  6     0.00000000E+00   # Real(ZE(3,6),dp)
  4  1    -9.99999997E-01   # Real(ZE(4,1),dp)
  4  2     1.25134556E-14   # Real(ZE(4,2),dp)
  4  3    -0.00000000E+00   # Real(ZE(4,3),dp)
  4  4     7.50775896E-05   # Real(ZE(4,4),dp)
  4  5     3.02207968E-18   # Real(ZE(4,5),dp)
  4  6    -0.00000000E+00   # Real(ZE(4,6),dp)
  5  1     1.25119018E-14   # Real(ZE(5,1),dp)
  5  2     9.99879579E-01   # Real(ZE(5,2),dp)
  5  3     0.00000000E+00   # Real(ZE(5,3),dp)
  5  4    -9.39133638E-19   # Real(ZE(5,4),dp)
  5  5    -1.55186263E-02   # Real(ZE(5,5),dp)
  5  6     0.00000000E+00   # Real(ZE(5,6),dp)
  6  1     0.00000000E+00   # Real(ZE(6,1),dp)
  6  2     0.00000000E+00   # Real(ZE(6,2),dp)
  6  3    -9.70485989E-01   # Real(ZE(6,3),dp)
  6  4     0.00000000E+00   # Real(ZE(6,4),dp)
  6  5     0.00000000E+00   # Real(ZE(6,5),dp)
  6  6     2.41157512E-01   # Real(ZE(6,6),dp)
Block SCALARMIX Q=  1.00000000E+03  # ()
  1  1    -6.38199818E-02   # ZH(1,1)
  1  2    -9.97961427E-01   # ZH(1,2)
  2  1    -9.97961427E-01   # ZH(2,1)
  2  2     6.38199818E-02   # ZH(2,2)
Block PSEUDOSCALARMIX Q=  1.00000000E+03  # ()
  1  1     6.23900398E-02   # ZA(1,1)
  1  2    -9.98051844E-01   # ZA(1,2)
  2  1    -9.98051844E-01   # ZA(2,1)
  2  2    -6.23900398E-02   # ZA(2,2)
Block CHARGEMIX Q=  1.00000000E+03  # ()
  1  1     6.23815999E-02   # Real(ZP(1,1),dp)
  1  2    -9.98052371E-01   # Real(ZP(1,2),dp)
  2  1    -9.98052371E-01   # Real(ZP(2,1),dp)
  2  2    -6.23815999E-02   # Real(ZP(2,2),dp)
Block NMIX Q=  1.00000000E+03  # ()
  1  1     9.99829929E-01   # Real(ZN(1,1),dp)
  1  2    -1.16851262E-03   # Real(ZN(1,2),dp)
  1  3     1.79404329E-02   # Real(ZN(1,3),dp)
  1  4    -4.10957336E-03   # Real(ZN(1,4),dp)
  2  1    -1.82117734E-03   # Real(ZN(2,1),dp)
  2  2    -9.99359235E-01   # Real(ZN(2,2),dp)
  2  3     3.36255202E-02   # Real(ZN(2,3),dp)
  2  4    -1.21295677E-02   # Real(ZN(2,4),dp)
  3  1     0.00000000E+00   # Real(ZN(3,1),dp)
  3  2    -0.00000000E+00   # Real(ZN(3,2),dp)
  3  3    -0.00000000E+00   # Real(ZN(3,3),dp)
  3  4    -0.00000000E+00   # Real(ZN(3,4),dp)
  4  1    -1.55428237E-02   # Real(ZN(4,1),dp)
  4  2     3.23762580E-02   # Real(ZN(4,2),dp)
  4  3     7.06387603E-01   # Real(ZN(4,3),dp)
  4  4    -7.06913540E-01   # Real(ZN(4,4),dp)
Block IMNMIX Q=  1.00000000E+03  # ()
  1  1     0.00000000E+00   # Aimag(ZN(1,1))
  1  2     0.00000000E+00   # Aimag(ZN(1,2))
  1  3     0.00000000E+00   # Aimag(ZN(1,3))
  1  4     0.00000000E+00   # Aimag(ZN(1,4))
  2  1     0.00000000E+00   # Aimag(ZN(2,1))
  2  2     0.00000000E+00   # Aimag(ZN(2,2))
  2  3     0.00000000E+00   # Aimag(ZN(2,3))
  2  4     0.00000000E+00   # Aimag(ZN(2,4))
  3  1     9.75792447E-03   # Aimag(ZN(3,1))
  3  2    -1.52161493E-02   # Aimag(ZN(3,2))
  3  3    -7.06798429E-01   # Aimag(ZN(3,3))
  3  4    -7.07184016E-01   # Aimag(ZN(3,4))
  4  1     0.00000000E+00   # Aimag(ZN(4,1))
  4  2     0.00000000E+00   # Aimag(ZN(4,2))
  4  3     0.00000000E+00   # Aimag(ZN(4,3))
  4  4     0.00000000E+00   # Aimag(ZN(4,4))
Block UMIX Q=  1.00000000E+03  # ()
  1  1     9.98870026E-01   # Real(UM(1,1),dp)
  1  2    -4.75254881E-02   # Real(UM(1,2),dp)
  2  1     4.75254881E-02   # Real(UM(2,1),dp)
  2  2     9.98870026E-01   # Real(UM(2,2),dp)
Block VMIX Q=  1.00000000E+03  # ()
  1  1     9.99852117E-01   # Real(UP(1,1),dp)
  1  2    -1.71972230E-02   # Real(UP(1,2),dp)
  2  1     1.71972230E-02   # Real(UP(2,1),dp)
  2  2     9.99852117E-01   # Real(UP(2,2),dp)
Block UELMIX Q=  1.00000000E+03  # ()
  1  1     1.00000000E+00   # Real(ZEL(1,1),dp)
  1  2     0.00000000E+00   # Real(ZEL(1,2),dp)
  1  3     0.00000000E+00   # Real(ZEL(1,3),dp)
  2  1     0.00000000E+00   # Real(ZEL(2,1),dp)
  2  2     1.00000000E+00   # Real(ZEL(2,2),dp)
  2  3     0.00000000E+00   # Real(ZEL(2,3),dp)
  3  1     0.00000000E+00   # Real(ZEL(3,1),dp)
  3  2     0.00000000E+00   # Real(ZEL(3,2),dp)
  3  3     1.00000000E+00   # Real(ZEL(3,3),dp)
Block UERMIX Q=  1.00000000E+03  # ()
  1  1     1.00000000E+00   # Real(ZER(1,1),dp)
  1  2     0.00000000E+00   # Real(ZER(1,2),dp)
  1  3     0.00000000E+00   # Real(ZER(1,3),dp)
  2  1     0.00000000E+00   # Real(ZER(2,1),dp)
  2  2     1.00000000E+00   # Real(ZER(2,2),dp)
  2  3     0.00000000E+00   # Real(ZER(2,3),dp)
  3  1     0.00000000E+00   # Real(ZER(3,1),dp)
  3  2     0.00000000E+00   # Real(ZER(3,2),dp)
  3  3     1.00000000E+00   # Real(ZER(3,3),dp)
Block UDLMIX Q=  1.00000000E+03  # ()
  1  1     1.00000000E+00   # Real(ZDL(1,1),dp)
  1  2     0.00000000E+00   # Real(ZDL(1,2),dp)
  1  3     0.00000000E+00   # Real(ZDL(1,3),dp)
  2  1     0.00000000E+00   # Real(ZDL(2,1),dp)
  2  2     1.00000000E+00   # Real(ZDL(2,2),dp)
  2  3     0.00000000E+00   # Real(ZDL(2,3),dp)
  3  1     0.00000000E+00   # Real(ZDL(3,1),dp)
  3  2     0.00000000E+00   # Real(ZDL(3,2),dp)
  3  3     1.00000000E+00   # Real(ZDL(3,3),dp)
Block UDRMIX Q=  1.00000000E+03  # ()
  1  1     1.00000000E+00   # Real(ZDR(1,1),dp)
  1  2     0.00000000E+00   # Real(ZDR(1,2),dp)
  1  3     0.00000000E+00   # Real(ZDR(1,3),dp)
  2  1     0.00000000E+00   # Real(ZDR(2,1),dp)
  2  2     1.00000000E+00   # Real(ZDR(2,2),dp)
  2  3     0.00000000E+00   # Real(ZDR(2,3),dp)
  3  1     0.00000000E+00   # Real(ZDR(3,1),dp)
  3  2     0.00000000E+00   # Real(ZDR(3,2),dp)
  3  3     1.00000000E+00   # Real(ZDR(3,3),dp)
Block UULMIX Q=  1.00000000E+03  # ()
  1  1     1.00000000E+00   # Real(ZUL(1,1),dp)
  1  2     0.00000000E+00   # Real(ZUL(1,2),dp)
  1  3     0.00000000E+00   # Real(ZUL(1,3),dp)
  2  1     0.00000000E+00   # Real(ZUL(2,1),dp)
  2  2     1.00000000E+00   # Real(ZUL(2,2),dp)
  2  3     0.00000000E+00   # Real(ZUL(2,3),dp)
  3  1     0.00000000E+00   # Real(ZUL(3,1),dp)
  3  2     0.00000000E+00   # Real(ZUL(3,2),dp)
  3  3     1.00000000E+00   # Real(ZUL(3,3),dp)
Block UURMIX Q=  1.00000000E+03  # ()
  1  1     1.00000000E+00   # Real(ZUR(1,1),dp)
  1  2     0.00000000E+00   # Real(ZUR(1,2),dp)
  1  3     0.00000000E+00   # Real(ZUR(1,3),dp)
  2  1     0.00000000E+00   # Real(ZUR(2,1),dp)
  2  2     1.00000000E+00   # Real(ZUR(2,2),dp)
  2  3     0.00000000E+00   # Real(ZUR(2,3),dp)
  3  1     0.00000000E+00   # Real(ZUR(3,1),dp)
  3  2     0.00000000E+00   # Real(ZUR(3,2),dp)
  3  3     1.00000000E+00   # Real(ZUR(3,3),dp)
Block SPhenoLowEnergy # low energy observables 
    1    3.26178118E-04  # BRBtoSGamma
   20    1.19535857E-14  # ae
   21    5.21991905E-10  # amu
   22    1.45303845E-07  # atau
   23    0.00000000E+00  # EDMe
   24    0.00000000E+00  # EDMmu
   25    0.00000000E+00  # EDMtau
   26    7.08322650E-33  # BRMuEgamma
   27    0.00000000E+00  # BRTauEgamma
   28    0.00000000E+00  # BRTauMugamma
   29    4.33984227E-35  # BrMu3e
   30    0.00000000E+00  # BrTau3e
   31    0.00000000E+00  # BrTau3mu
   39    5.69929370E-05  # dRho
   40    1.64302946E-45  # BRZMuE
   41    0.00000000E+00  # BRZTauE
   42    0.00000000E+00  # BRZTauMu
   43    8.44021768E-14  # BRBsEE
   44    3.61252520E-09  # BRBsMuMu
   45    0.00000000E+00  # BRBsMuE
   46    2.77043277E-15  # BRBdEE
   47    1.18573162E-10  # BRBdMuMu
   48    2.48465147E-08  # BRBdTauTau
   81    1.93168683E-35  # MuEAl
   82    3.47676949E-35  # MuETi
   83    4.70560213E-35  # MuESr
   84    5.29337489E-35  # MuESb
   85    2.85151309E-35  # MuEAu
   86    2.68223255E-35  # MuEPb
   91    0.00000000E+00  # TauEPi0
   92    0.00000000E+00  # TauEEta
   93    0.00000000E+00  # TauEEtap
   94    0.00000000E+00  # TauMuPi0
   95    0.00000000E+00  # TauMuEta
   96    0.00000000E+00  # TauMuEtap
  800    2.35135983E+03  # psMuEAl
  801    4.71491063E+03  # psMuETi
  802    5.90097458E+03  # psMuESr
  803    8.61708893E+03  # psMuESb
  804    5.17097974E+03  # psMuEAu
  805    4.99448558E+03  # psMuEPb
Block SPheno # SPheno internal parameters 
         1    -1.00000000E+00  # ErrorLevel
         2     1.00000000E+00  # SPA_conventions
        11     0.00000000E+00  # Branching ratios
        13     1.00000000E+00  # 3 Body decays
        31     1.34632464E+16  # GUT scale
        33     1.00000000E+03  # SUSY scale
        34     1.00000000E-04  # Precision
        35     4.00000000E+01  # Iterations
        38     2.00000000E+00  # RGE level
        40     7.29735308E-03  # Alpha
        41     2.49000000E+00  # Gamma_Z
        42     2.06000000E+00  # Gamma_W
        50     1.00000000E+00  # Rotate negative fermion masses
        51     0.00000000E+00  # Switch to SCKM matrix
        52     0.00000000E+00  # Ignore negative masses
        53     0.00000000E+00  # Ignore negative masses at MZ
        55     1.00000000E+00  # Calculate one loop masses
        56     0.00000000E+00  # Calculate two-loop Higgs masses
        57     1.00000000E+00  # Calculate low energy
        60     1.00000000E+00  # Include kinetic mixing
        65     1.00000000E+00  # Solution of tadpole equation
BLOCK VEVACIOUSRESULTS # results from Vevacious version 0.3.0, documented in [still unpublished]
    0   0    -1.00000000E+000    short-lived    # stability of input
    0   1    +7.24581244E-020    direct_path_bound    # tunneling time in Universe ages / calculation type
    1   0    -1.09640888E+008    relative_depth    # input potential energy density difference from all VEVs = 0.0
    1   1    +0.00000000E+000    vE3    # input VEV
    1   2    +0.00000000E+000    vL3    # input VEV
    1   3    +2.27373675E-010    vQ3    # input VEV
    1   4    +0.00000000E+000    vT3    # input VEV
    1   5    +1.52965823E+001    vdR    # input VEV
    1   6    +2.43954795E+002    vuR    # input VEV
    2   0    -3.80042153E+014    relative_depth    # global minimum potential energy density difference from all VEVs = 0.0
    2   1    +1.42993056E+004    vE3    # global minimum VEV
    2   2    +1.41690972E+004    vL3    # global minimum VEV
    2   3    +5.09468337E+003    vQ3    # global minimum VEV
    2   4    +5.19310255E+003    vT3    # global minimum VEV
    2   5    +1.47001915E+004    vdR    # global minimum VEV
    2   6    +5.78713893E+003    vuR    # global minimum VEV
