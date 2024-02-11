function [f,ceq]=mycon1(x,btk,stk,snum,bnum,stpv,btpv,Esd)

f=[btk/x(1)+stk/x(2)-snum-bnum-stpv-btpv+Esd];
ceq=[];

end