function apa()
	%generera tridiagonala matrisen
	A = triDiag(100);
	
	% homogenledet
	noll = zeros(100,1);
	% hastighetsledet
	v = zeros(100,1);
	%ger mitterta partiklarna en hastighet vid t=0
	v(45:55,1) = ones(11,1);
	v = 0.1 * v
	
	rhs(1:100, 1) = noll;
	rhs(1:100, 2) = v;
	
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
	
	[P D]  = eig(A);
	lambda = D * ones(100, 1);
	
	%skriver högerleden som en konstanter
	solution = P \ rhs;
	
	a = solution(:, 1);
	b = solution(:, 2);
	
	%vi låter omega0 vara tidsenhet, (fås hastighet i enhet meter*omega0)
	k = sqrt(lambda); %*omega0
	C = b./k;
	%C = sqrt(c.^2+a.^2);
	fi = 0;%fi = asin(a./C); %blir 0 alltid, ty Aa=0 saknar trivial lösning, enl. IMT; dock: ±n*pi
	
	%ger förenklad lösning: C=c=k./b=sqrt(lambda)./(A\v)
	result = @(t) [P * (C .* sin(sqrt(lambda) .* t + fi))];
	
	time  = 150;
	stime = 0;
	
	num_steps = 200;
	time_step = time / num_steps;
	
	for t=0:num_steps-1
		data(1:100, t+1) = result(t*time_step + stime);
		
		h = plot(linspace(0, 100), result(t*time_step + stime));
		ylim([-1 1]);
		saveas(h, strcat('plot', sprintf('%d', t), '.png'));
	end
	
	plot(repmat(linspace(stime, stime+time, num_steps), 100, 1)', data)
	
	%t = linspace(0,10); %kommer behöva ändras för andra storlekar på A, 100 element just nu
	%tplot = repmat(t,100,1); %repmat för dimension
	%result(t')
	%plot(tplot',result(t')); %t transponat för behöver kolonn
	%plotten ska tolkas hur? ser konstigt ut?
function [matrix] = triDiag(side_length)
	%Generera den tridiagonala matrisen:
	n = -ones(side_length - 1, 1);
	B = diag(n, 1);
	C = diag(n, -1);
	n = 2 * ones(side_length, 1);
	A = diag(n);
	matrix = A + B + C;