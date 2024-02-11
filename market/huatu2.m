k=50;
pv=20;
pb=0.4;
ps=0.8;
a=0.42;
x=0:0.1:200;


for i=1:2001
    pro(i)=k*log(1+x(i))-ps*max(x(i)-pv,0)-pb*min(x(i)-pv,0)+pv*a;
end
[z,h]=max(pro);
plot(0:0.1:200,pro);
hold on
plot(h*0.1,z,'*');

hold on
ps=0.9;
for i=1:2001
    pro(i)=k*log(1+x(i))-ps*max(x(i)-pv,0)-pb*min(x(i)-pv,0)+pv*a;
end

[z,h]=max(pro);
plot(0:0.1:200,pro);
hold on
plot(h*0.1,z,'*');

hold on
ps=1;
for i=1:2001
    pro(i)=k*log(1+x(i))-ps*max(x(i)-pv,0)-pb*min(x(i)-pv,0)+pv*a;
end

[z,h]=max(pro);
plot(0:0.1:200,pro);
hold on
plot(h*0.1,z,'*');

hold on
k=30;
pv=80;
pb=0.4;
ps=1;
for i=1:2001
    pro(i)=k*log(1+x(i))-ps*max(x(i)-pv,0)-pb*min(x(i)-pv,0)+pv*a;
end

[z,h]=max(pro);
plot(0:0.1:200,pro);
hold on
plot(h*0.1,z,'*');

hold on
pb=0.5;
for i=1:2001
    pro(i)=k*log(1+x(i))-ps*max(x(i)-pv,0)-pb*min(x(i)-pv,0)+pv*a;
end

[z,h]=max(pro);
plot(0:0.1:200,pro);
hold on
plot(h*0.1,z,'*');

hold on
pb=0.6;
for i=1:2001
    pro(i)=k*log(1+x(i))-ps*max(x(i)-pv,0)-pb*min(x(i)-pv,0)+pv*a;
end

[z,h]=max(pro);
plot(0:0.1:200,pro);
hold on
plot(h*0.1,z,'*');