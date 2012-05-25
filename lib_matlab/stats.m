function[] = stats(n1)

IN= '../CHAnWillB/OUT/stats.dat'

figure(1);
clf;
figure(2);
clf;
theend = logical(0);
f = fopen(IN,'r');
hold on;

[A]=fscanf(f,'%f', [3,n1]);

fclose(f);

for i=1:n1
  ta(i) = A(1,i);
  ma(i) = A(2,i);
  ea(i) = A(3,i);
end;

maxmass = max(ma)
minmass = min(ma)
trumass = ma(1)

massdiff = (maxmass-minmass)/trumass

figure(1);
plot(ta,ma)

hold on;

figure(2)
plot(ta,ea)

hold on;


