function [f,ceq]=mycon1(x,btk,stk,kall,snum,bnum,stpv,btpv,Esd)

f=[sum(kall(2,:)./x(1,:))+sum(kall(3,:)./x(2,:))-snum-bnum-stpv-btpv+Esd];
ceq=[];

end