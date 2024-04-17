clc
clear

%input section
%inputs mesh
xnodes = linspace(0,8,33);
GL = 4;
n = size(xnodes,2)-1;
for e = 1:n
    L(e) = xnodes(e+1)-xnodes(e);
end

syms u [2*(n+1),1];
syms r [2*(n+1),1];
syms x;

%input boundary displacements
u(1) = 0;
u(33) = 0;
u(65) = 0;

%moment of inertia
for e = 1:n
    EI(e) = 2E4;
end

for i = 1:(2*(n+1))
    if u(i) ~= 0
        r(i) = 0;
    end
end

%input boundry and field forces
P = [2 -10];
M = [6 5];
px = [0 8 -2];

%end of input section

%gather matrix
R = zeros(GL,2*(n+1),n);
for e = 1:n
    for i = 1:GL
        for j = 1:(2*n+2)
            if j == i+2*(e-1)
            R(i,j,e) = 1;
            end
        end
    end
end
Rt = pagetranspose(R);

%shape functions
for e = 1:n
    Nu1(x) = 1/4*(1-(2*(x-xnodes(e))/L(e)-1))^2*(2+(2*(x-xnodes(e))/L(e)-1));
    Nt1(x) = L(e)/8*(1-(2*(x-xnodes(e))/L(e)-1))^2*(1+(2*(x-xnodes(e))/L(e)-1));
    Nu2(x) = 1/4*(1+(2*(x-xnodes(e))/L(e)-1))^2*(2-(2*(x-xnodes(e))/L(e)-1));
    Nt2(x) = L(e)/8*(1+(2*(x-xnodes(e))/L(e)-1))^2*((2*(x-xnodes(e))/L(e)-1)-1);
    Ne(:,:,e) = [Nu1 Nt1 Nu2 Nt2];
end

%stiffness matrix
Ke = zeros(GL,GL,n);
for e = 1:n
    Ke(:,:,e) = int(diff(diff(Ne(:,:,e)'))*EI(e)*diff(diff(Ne(:,:,e))),x,xnodes(e),xnodes(e+1));
end
K = zeros(2*(n+1));
for e = 1:n
    K = K + Rt(:,:,e)*Ke(:,:,e)*R(:,:,e);
end

%element boundary force matrix
for e = 1:n
    fGm(:,:) = zeros(GL,n);
    fGs(:,:) = zeros(GL,n);
end
for e = 1:n
    N(x) = Ne(:,:,e);
    for i = 1:size(P,1)

        if P(i,1) == xnodes(e) || (P(i,1) == xnodes(n+1) && e == n)
            fGs(:,e) = N(P(i,1))'*P(i,2);

        elseif P(i,1) > xnodes(e) && P(i,1) < xnodes(e+1)
            fGs(:,e) = N(P(i,1))'*P(i,2);

        end
    end
    N(x) = diff(Ne(:,:,e),x);
    for i = 1:size(M,1)

        if M(i,1) == xnodes(e) || (M(i,1) == xnodes(n+1) && e == n)
            fGm(:,e) = N(M(i,1))'*M(i,2);

        elseif M(i,1) > xnodes(e) && M(i,1) < xnodes(e+1)
            fGm(:,e) = N(M(i,1))'*M(i,2);
            
        end
    end
end

%element field force matrix
for e = 1:n
fO(:,:) = zeros(GL,n);
end
for e = 1:n
    for i = 1:size(px,1)

        if px(i,1) > xnodes(e) && px(i,2) < xnodes(e+1)
            fO(:,e) = int((Ne(:,:,e)')*px(i,3),x,px(i,1),px(i,2));

        elseif px(i,1) <= xnodes(e) && px(i,2) >= xnodes(e+1)
            fO(:,e) = int((Ne(:,:,e)')*px(i,3),x,xnodes(e),xnodes(e+1));

        elseif (px(i,1) > xnodes(e) && px(i,1) < xnodes(e+1)) && (px(i,2) >= xnodes(e+1))
            fO(:,e) = int((Ne(:,:,e)')*px(i,3),x,px(i,1),xnodes(e+1));

        elseif px(i,1) <= xnodes(e) && (px(i,2) < xnodes(e+1) && px(i,2) > xnodes(e))
            fO(:,e) = int((Ne(:,:,e)')*px(i,3),x,xnodes(e),px(i,2));

        end
    end
end

%global force matrix
f = zeros(2*(n+1),1);
for e = 1:n
    f = f + Rt(:,:,e)*fGs(:,e) + Rt(:,:,e)*fGm(:,e) + Rt(:,:,e)*fO(:,e);
end

%solver
j = 1;
for i = 1:2*(n+1)
    if symType(u(i))~="integer"
        solve_u(j) = u(i);
        j = j + 1;
    end
end
j = 1;
for i = 1:2*(n+1)
    if symType(r(i))~="integer"
        solve_r(j) = r(i);
        j = j + 1;
    end
end

seq = cat(2,solve_u,solve_r);
solution = vpasolve(K*u==f+r,seq)

fields = fieldnames(solution);

for i = 1:numel(fields)
    s(i) = solution.(fields{i});
end

k = 1;
for i = 1:2*(n+1)
    if symType(u(i))~="integer"
        u(i) = s(k);
        k = k + 1;
    end
end

d = double(u)




%draw graphs
step = 0.001;

RA = 447/64;
RB = 585/32;
RC = 47/64;
C = -239/24;

% exact displacement field
x = 0:step:xnodes(n+1);
yex =((-(x.^4)/12 + (RA/6)*x.^3 - (10/6)*(mac(x-2).^3) + (RB/6)*(mac(x-4).^3) - (5/2)*(mac(x-6).^2) + (C*x))/EI(e))';

% fe displacement field
yfe = [0];
for e = 1:n
    xe = 0:step:L(e);
    Nu1 = (((2*xe)/L(e) + 1).*((2*xe)/L(e) - 2).^2)/4;
    Nt1 = (xe.*((2*xe)/L(e) - 2).^2)/4;
    Nu2 = -(xe.^2.*((2*xe)/L(e) - 3))/L(e)^2;
    Nt2 = (xe.^2.*((2*xe)/L(e) - 2))/(2*L(e));
    Ne = [Nu1' Nt1' Nu2' Nt2'];
    Ned = Ne*R(:,:,e);
    plt = Ned*d;
    yfe(size(yfe,1),:) = [];
    yfe = cat(1,yfe,plt);
end

%exact moment field
x = 0:step:xnodes(n+1);
mex = ((-(x.^2) + (RA)*x - (10)*(mac(x-2)) + (RB)*(mac(x-4)) - (5)*(macm(x-6))))';

%fe moment field
mfe = [0];
for e = 1:n
    xe = 0:step:L(e);
    Nu1 = (2*((2*xe)/L(e) + 1))/L(e)^2 + (4*((2*xe)/L(e) - 2))/L(e)^2;
    Nt1 = (2*xe)/L(e)^2 + (2*((2*xe)/L(e) - 2))/L(e);
    Nu2 = - (8*xe)/L(e)^3 - (2*((2*xe)/L(e) - 3))/L(e)^2;
    Nt2 = (4*xe)/L(e)^2 + ((2*xe)/L(e) - 2)/L(e);
    Ne = [Nu1' Nt1' Nu2' Nt2'];
    Ned = Ne*R(:,:,e);
    plt = EI(e)*Ned*d;
    mfe(size(mfe,1),:) = [];
    mfe = cat(1,mfe,plt);
end

%exact shear field
x = 0:step:xnodes(n+1);
vex = ((-(2*x) + (RA) - (10)*(macm(x-2)) + (RB)*(macm(x-4))))';

%fe shear field
vfe = [0];
for e = 1:n
    xe = 0:step:L(e);
    Nu1 = 12/L(e)^3*(ones(1,size(xe,2)));
    Nt1 = 6/L(e)^2*(ones(1,size(xe,2)));
    Nu2 = -12/L(e)^3*(ones(1,size(xe,2)));
    Nt2 = 6/L(e)^2*(ones(1,size(xe,2)));
    Ne = [Nu1' Nt1' Nu2' Nt2'];
    Ned = Ne*R(:,:,e);
    plt = EI(e)*Ned*d;
    vfe(size(vfe,1),:) = [];
    vfe = cat(1,vfe,plt);
end

fig = figure(1)

%draw displacements
subplot(3,1,1)
plot(x,1000*yex,'black')
hold on
plot(x,1000*yfe,'red')
xlabel('x (m)')
ylabel('deflexão (mm)')
title(sprintf('Deflexão com %d elementos',n))

%draw moments
subplot(3,1,2)
plot(x,mex,'black')
hold on
plot(x,mfe,'red')
xlabel('x (m)')
ylabel('momento fletor (kNm)')
title(sprintf('Momento fletor com %d elementos',n))

%draw shear force
subplot(3,1,3)
plot(x,vex,'black')
hold on
plot(x,vfe,'red')
hold on
xlabel('x (m)')
ylabel('força cortante (kN)')
title(sprintf('Força cortante com %d elementos',n))

legend('Solução Exata','Aproximação MEF','Location','southeast')

function res = mac(ch)
    for i = 1:size(ch,2)
        if ch(i) <=0
            res(i) = 0;
        else
            res(i) = ch(i);
        end
    end
end

function res = macm(ch)
    for i = 1:size(ch,2)
        if ch(i) <=0
            res(i) = 0;
        else
            res(i) = 1;
        end
    end
end
