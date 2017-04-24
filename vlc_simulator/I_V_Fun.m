function IF = I_V_Fun(VF,VT,n,ISat)









aux = VF > 0;




IF = ISat.*(exp(VF./(n*VT)) - 1) .* aux;
