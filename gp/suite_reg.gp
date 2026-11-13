default(realprecision, 60)
foreach([41,65,137,145,161,219,257,265,291,299,313,323,353,371,395,457,465,505,561,609,651,689,721,723,731,761,777,793,889,905,915,985,987,995], d, E = ellinit([0,0,0,-d^2,0]); rk = ellrank(E); P = rk[4][1]; Q = rk[4][2]; hp = ellheight(E,P); hq = ellheight(E,Q); s = ellheight(E,elladd(E,P,Q))-hp-hq; reg = hp*hq-s^2/4; print("G ",d,";",rk[1],";",rk[2],";",reg,";",lfun(E,1,2)))
quit()
