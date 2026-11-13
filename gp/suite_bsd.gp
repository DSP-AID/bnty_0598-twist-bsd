D = 1000
h = 10.^-13
forstep(d=1, D, 2, if(issquarefree(d), E = ellinit([0,0,0,-d^2,0]); w = ellrootno(E); V = if(w == 1, lfun(E, 1), (8*(lfun(E,1+h)-lfun(E,1-h))-(lfun(E,1+2*h)-lfun(E,1-2*h)))/(12*h)); print("B ", d, " ", ellglobalred(E)[1], " ", Str(ellperiods(E)[1]), " ", elltamagawa(E), " ", Str(elltors(E)[2]), " ", w, " ", V)))
quit()
