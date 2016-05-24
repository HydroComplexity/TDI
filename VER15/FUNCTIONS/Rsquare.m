%---------The function return R-square value -------------
%--------- input data ==> x-data, y-data and p-data ----
%--------- output data ==> r2
function [r2]=Rsquare(x,y,p)

Ymeasure=y;
Ycalculate=(p(1)*x)+p(2);
meanY=mean(Ymeasure);
deltaY2=sum((Ycalculate-Ymeasure).^2);
distanceY2=sum((Ymeasure-meanY).^2);
r2=1-(deltaY2/distanceY2);
end
%-------------------------------------------------------