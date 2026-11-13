D = 1000
forstep(d=1, D, 2, if(issquarefree(d), E = ellinit([0,0,0,-d^2,0]); print("R ", d, " ", ellrootno(E))))
quit()
