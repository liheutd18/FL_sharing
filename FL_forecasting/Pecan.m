function [Newntl,Newl,Newpv]=Pecan()

    load pecanYL.mat
    for i=1:20
        NTL(:,i)=PecanYL(:,2*i-1);
        PV(:,i)=PecanYL(:,2*i);
    end
    NTL=NTL(:,5:15);
    PV=PV(:,5:15);
    PV(338,1)=3.40;
    PV(find( PV<0))=0;
    %estimatepv(10,:)=0;estimatepv(8,:)=0;estimatepv(5,:)=0;

    NTL(:,6)=[];PV(:,6)=[];
    PV(:,10)=0;PV(:,10)=0;PV(:,10)=0;
    

    L=NTL+PV;
    L=abs(L);
    k=length(NTL);
 
    for i=1:10
         a=reshape(NTL(:,i),4,k/4);
         b=reshape(PV(:,i),4,k/4);
         c=reshape(L(:,i),4,k/4);
         Newntl(:,i)=sum(a)'/4; 
         Newpv(:,i)=sum(b)'/4; 
         Newl(:,i)=sum(c)'/4; 
    end


    pv=Newpv;
    l=Newl;
    
    kmin=0.8;
    kmax=1.2;
    loadmin=kmin*l;
    loadmax=kmax*l;

end
    
