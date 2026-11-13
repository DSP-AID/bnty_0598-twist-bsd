D = 250
forstep(d=1, D, 2, if(issquarefree(d), E = ellinit([0,0,0,-d^2,0]); rk = ellrank(E); if(rk[1] > 0, print("K ", d, " ", rk[1], " ", Str([ellheight(E, P) | P <- rk[4]])), print("K ", d, " 0"))))
quit()
