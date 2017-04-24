function G = freqRespLED(w)

wc = 2*pi*1e6;

w2 = 2*pi*3.26e6;

w3 = 2*pi*10.86e6;


wAux1 = w(abs(w)<=wc);

wAux2 = w(abs(w)>wc);

G = zeros(length(wAux1) + length(wAux2),1);

G(abs(w)<=wc) = exp(-abs(wAux1)/w2);

G(abs(w)>wc) = exp(-abs(wc)/w2)*exp(abs(wc)/w3)*exp(-abs(wAux2)/w3);

