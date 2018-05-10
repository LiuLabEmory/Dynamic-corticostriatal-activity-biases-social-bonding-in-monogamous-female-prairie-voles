function Sc = lnrusa(S,Fs)

bns = length(S);
tv = (0:bns-1)'/Fs;
Sc = S;
frqlnr = [59.8;59.9;60.0;60.1;60.2;119.8;119.9;120.0;120.1;120.2;179.8;179.9;180.0;180.1;180.2];
%frqlnr = 59:0.1:61;

for lop = 1:length(frqlnr)
    frqlop = frqlnr(lop);
    Sft = 2*sum(S .* exp(i*2*pi*frqlop*tv))/bns;
    Sc = Sc - abs(Sft)*cos(2*pi*frqlop*tv - angle(Sft));
end   