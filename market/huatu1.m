k=10;
pv=40;
ps=0.9;
pb=0.5;
a=0.42;
x=0:0.1:200;
h=zeros(1,10);
z=zeros(1,10);

for i=1:2001
    pro(i)=k*log(1+x(i))-ps*max(x(i)-pv,0)-pb*min(x(i)-pv,0)+pv*a;
end
[z,h]=max(pro);
plot(0:0.1:200,pro);
hold on
plot(h*0.1,z,'*');

hold on
k=20;
for i=1:2001
    pro(i)=k*log(1+x(i))-ps*max(x(i)-pv,0)-pb*min(x(i)-pv,0)+pv*a;
end

[z,h]=max(pro);
plot(0:0.1:200,pro);
hold on
plot(h*0.1,z,'*');

hold on
k=30;
for i=1:2001
    pro(i)=k*log(1+x(i))-ps*max(x(i)-pv,0)-pb*min(x(i)-pv,0)+pv*a;
end

[z,h]=max(pro);
plot(0:0.1:200,pro);
hold on
plot(h*0.1,z,'*');

hold on
k=40;
for i=1:2001
    pro(i)=k*log(1+x(i))-ps*max(x(i)-pv,0)-pb*min(x(i)-pv,0)+pv*a;
end

[z,h]=max(pro);
plot(0:0.1:200,pro);
hold on
plot(h*0.1,z,'*');

hold on
k=50;
for i=1:2001
    pro(i)=k*log(1+x(i))-ps*max(x(i)-pv,0)-pb*min(x(i)-pv,0)+pv*a;
end

[z,h]=max(pro);
plot(0:0.1:200,pro);
hold on
plot(h*0.1,z,'*');

hold on
k=60;
for i=1:2001
    pro(i)=k*log(1+x(i))-ps*max(x(i)-pv,0)-pb*min(x(i)-pv,0)+pv*a;
end

[z,h]=max(pro);
plot(0:0.1:200,pro);
hold on
plot(h*0.1,z,'*');

hold on
k=70;
for i=1:2001
    pro(i)=k*log(1+x(i))-ps*max(x(i)-pv,0)-pb*min(x(i)-pv,0)+pv*a;
end

[z,h]=max(pro);
plot(0:0.1:200,pro);
hold on
plot(h*0.1,z,'*');

hold on
k=80;
for i=1:2001
    pro(i)=k*log(1+x(i))-ps*max(x(i)-pv,0)-pb*min(x(i)-pv,0)+pv*a;
end

[z,h]=max(pro);
plot(0:0.1:200,pro);
hold on
plot(h*0.1,z,'*');

hold on
k=90;
for i=1:2001
    pro(i)=k*log(1+x(i))-ps*max(x(i)-pv,0)-pb*min(x(i)-pv,0)+pv*a;
end

[z,h]=max(pro);
plot(0:0.1:200,pro);
hold on
plot(h*0.1,z,'*');

hold on
k=100;
for i=1:2001
    pro(i)=k*log(1+x(i))-ps*max(x(i)-pv,0)-pb*min(x(i)-pv,0)+pv*a;
end

[z,h]=max(pro);
plot(0:0.1:200,pro);
hold on
plot(h*0.1,z,'*');