function[] = meshcntrplot2d(s0,nn,var,toggle0,toggle1,toggle2)
%
% s0 is the string indicating the appropriate PROBLEM directory, eg, 'CH'
% nn is the frame number to be printed
% var is the variable number
% toggle0 = 0 for adaptive mesh
% toggle0 = 1 for uniform mesh
% toggle1 = 0 for no output
% toggle1 = 1 for jpg output
% toggle1 = 2 for eps output
% toggle2 = 0 for plot only on the fine level patches, showing the mesh
% toggle2 = 1 for use pcolor on all level patches without showing the mesh

s1 = ['0000000' num2str(nn)];
s2 = s1((length(s1)-4):length(s1));

dir =['../' s0 '/OUT/']

if toggle0 == 0
  s3 = 'm'
else
  s3 = 'u'
end;

IN  = [dir s3 s2 '.dat']
if toggle1 == 1
  OUT = [dir s3 s2 '.jpg']
elseif toggle1 == 2
  OUT = [dir s3 s2 '.eps']
end;

clf reset;

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.5 0.5 4 4]);

theend = logical(0);
f = fopen(IN,'r');
hold on;
ipatch=0;

[time,count] = fscanf(f, '%f', 1);
[maxlevel,count] = fscanf(f, '%d', 1);
   
while(~theend)

  [level,count] =  fscanf(f, '%d', 1);
  [ndim,count]  =  fscanf(f, '%d', 1);
  [r,count] =  fscanf(f, '%d', 1);
  [nrvars,count]  =  fscanf(f, '%d', 1);

  if count ~= 0

    [dx(1),count] =  fscanf(f, '%f', 1);
    [dx(2),count] =  fscanf(f, '%f', 1);

    [xl(1),count] =  fscanf(f, '%f', 1);
    [xl(2),count] =  fscanf(f, '%f', 1);

    [xu(1),count] =  fscanf(f, '%f', 1);
    [xu(2),count] =  fscanf(f, '%f', 1);

    [n(1),count] =  fscanf(f, '%d', 1);
    [n(2),count] =  fscanf(f, '%d', 1);
    
    [mg(1,1),count] =  fscanf(f, '%d', 1);
    [mg(1,2),count] =  fscanf(f, '%d', 1);
    [mg(2,1),count] =  fscanf(f, '%d', 1);
    [mg(2,2),count] =  fscanf(f, '%d', 1);

    ipatch = ipatch + 1;
    disp(sprintf('Processing patch number %d', ipatch));
    disp(sprintf('The patch size is %d,  %d', n(1), n(2)));

    xu = xl+dx.*n;
    if toggle2 == 0
      drawmeshpatch(n,xl,xu,dx);
    end;
    
    if level == 0
      xlg(1) = xl(1);
      xug(1) = xu(1);
      xlg(2) = xl(2);
      xug(2) = xu(2);
    end;

    A = zeros(nrvars,(n(1)+2)*(n(2)+2),'double');
    
    [A]=fscanf(f,'%f', [nrvars,(n(1)+2)*(n(2)+2)]); % ghost layer included.
   
    if or(level == maxlevel,toggle2 == 1)
      drawcont(A,var,n,xl,xu,dx,xlg,xug,toggle2);
    end;
  else
    theend = 1;
  end
end;

text(xug(1)-1.4,xug(2)+0.2,['time = ' num2str(time)], 'FontSize', 14);

grid on
t = (xug-xlg)./4;
set(gca,'XTick',xlg(1):t(1):xug(1))
set(gca,'YTick',xlg(2):t(2):xug(2))

if toggle1 == 1
  print('-djpeg','-r400',OUT)
elseif toggle1 == 2
  print('-deps','-r1200',OUT)
end;

ipatch
fclose(f);

function drawcont(A,var,n,xl,xu,dx,xlg,xug,toggle2)

c = zeros(n(1)+2,n(2)+2,'double');
x = zeros(n(1)+2,n(2)+2,'double');
y = zeros(n(1)+2,n(2)+2,'double');
x1 = zeros(n(1)/2+1,n(2)/2+1,'double');
y1 = zeros(n(1)/2+1,n(2)/2+1,'double');
u = zeros(n(1)/2+1,n(2)/2+1,'double');
v = zeros(n(1)/2+1,n(2)/2+1,'double');

xll = xl(1)-dx(1)/2.0;
yll = xl(2)-dx(2)/2.0;

for j = 1:n(2)+2
  for i = 1:n(1)+2
    x(i,j) = xll+(i-1)*dx(1);
    y(i,j) = yll+(j-1)*dx(2);
    c(i,j) = A(var,(j-1)*(n(1)+2)+i);
    p(i,j) = A(4,(j-1)*(n(1)+2)+i);
  end;
end;

for j = 1:n(2)/2+1
  for i = 1:n(1)/2+1
    x1(i,j) = (x(2*i,2*j)+x(2*i-1,2*j)+x(2*i,2*j-1)+x(2*i-1,2*j-1))/4;
    y1(i,j) = (y(2*i,2*j)+y(2*i-1,2*j)+y(2*i,2*j-1)+y(2*i-1,2*j-1))/4;
    u(i,j) = -(p(2*i,2*j)-p(2*i-1,2*j)+p(2*i,2*j-1)-p(2*i-1,2*j-1))/4;
    v(i,j) = -(p(2*i,2*j)+p(2*i-1,2*j)-p(2*i,2*j-1)-p(2*i-1,2*j-1))/4;
  end;
end;

if toggle2 == 0

  %contour(x,y,c,0.5)
  %colormap([0,0,0])
  %contour(x,y,c,0.5,'Linewidth',1.5)
  %colormap([0,0,0])
  %colormap('default')
  w = [.5];
  %w = [.1 .3 .5 .7 .9];
  %contourf(x,y,c,w,'LineStyle','none')
  %contourf(x,y,c)
  hold on
  contour(x,y,c,w,'Linewidth',1.5)
  hold on
  quiver(x1,y1,u,v,2)
  hold on
  %contour(x,y,c,v,'Color','black','LineWidth',1)
  %[C,h] = contour(x,y,c,20);
  %clabel(C,h);
  %colormap('default')
  %surf(x,y,c,'LineStyle','none')

elseif toggle2 ==1

  pcolor(x,y,c)
  shading flat
  colormap(bone)
  colorbar;

end;

%pause

axis([xlg(1),xug(1),xlg(2),xug(2)])
axis equal
axis([xlg(1),xug(1),xlg(2),xug(2)])

function[] = drawmeshpatch(n,xl,xu,dx)

line([xl(1),xu(1)],[xl(2),xl(2)],'Color','black','LineWidth',1.5)
line([xl(1),xu(1)],[xu(2),xu(2)],'Color','black','LineWidth',1.5)

for j=1:n(2)-1
  y = xl(2)+j*dx(2);
  line([xl(1),xu(1)],[y,y],'Color','black','LineWidth',1.0)
end;

line([xl(1),xl(1)],[xl(2),xu(2)],'Color','black','LineWidth',1.5)
line([xu(1),xu(1)],[xl(2),xu(2)],'Color','black','LineWidth',1.5)

for i=1:n(1)-1
  x = xl(1)+i*dx(1);
  line([x,x],[xl(2),xu(2)],'Color','black','LineWidth',1.0)
end;
