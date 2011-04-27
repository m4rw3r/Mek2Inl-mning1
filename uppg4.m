function apa()
%generera tridiagonala matrisen
A=triDiag(100);
%homogenledet
noll=zeros(100,1);
%hastighetsledet
v=zeros(100,1);
%ger mitterta partiklarna en hastighet vid t=0
v(45:55,1)=ones(11,1);

%skriver högerleden som en konstanter
a=A\noll;
b=A\v;
%__________
%Vi har ekvationerna:
%A C.*sin(fi)=noll                      (1)
%A sqrt(lambda)*omega0.*C.*cos(Fi)=v    (2)
%(1) ger:
%=>  fi=arcsin((A\noll)/C)=arcsin(a/C)
%<=> arccos(sqrt(C.^2+a.^2)/C)
%(2) ger:
%A*C.*cos(Fi)=v./(sqrt(lambda)*omega0)
%C.*cos(Fi)=A\(v./(sqrt(lambda)*omega0))
%Varpå C fås från insättning av Fi, från (1), och Fi fås från C.
%___________
[P D]=eig(A);
lambda=D*ones(100,1);
%vi låter omega0 vara tidsenhet, (fås hastighet i enhet meter*omega0)
k=sqrt(lambda); %*omega0
c=k./b;
C=sqrt(c.^2+a.^2);
fi=asin(a./C); %blir 0 alltid, Aa=0 saknar trivial lösning? även ±n*pi? lurigt
result = @(t) [A * (C .* sin(sqrt(lambda).*t + fi))];
bla=A * (C .* sin(sqrt(lambda) + fi))

t=linspace(0,10); %kommer behöva ändras för andra storlekar på A, 100 element just nu
tplot=repmat(t,100,1); %repmat för dimension
plot(tplot,result(t')); %transponat för behöver kolonn

function [matrix] = triDiag(side_length)
	%Generera den tridiagonala matrisen:
	n = -ones(side_length - 1, 1);
	B = diag(n, 1);
	C = diag(n, -1);
	n = 2 * ones(side_length, 1);
	A = diag(n);
	matrix = A + B + C;