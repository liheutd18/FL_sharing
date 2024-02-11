function f=myfun2(x,stk,btk,snum,bnum,stpv,btpv,pgb,Esd)

f=stk+btk-x(1)*(bnum+btpv)-x(2)*(snum+stpv)+pgb*(stpv+btpv+snum+bnum-(btk/x(1))-(stk/x(2))-Esd);
f=-f;

end