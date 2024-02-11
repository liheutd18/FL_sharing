function f=myfun3future(price,stk,btk,snum,bnum,stpv,btpv,pgs,Esd)
ab=length(stk);
x1=price(1)*ones(ab,1);x2=price(2)*ones(ab,1);
f=stk+btk-x1.*(bnum+btpv)-x2.*(snum+stpv)+pgs*(stpv+btpv+snum+bnum-(btk./x1)-(stk./x2)-Esd);
f=sum(-f);
end