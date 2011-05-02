function apa()
clf
figure(1)
clf
	%generera tridiagonala matrisen
	A = triDiag(100);
	
	% homogenledet
	noll = zeros(100,1);
	% hastighetsledet
	v = zeros(100,1);
	%ger mitterta partiklarna en hastighet vid t=0
	v(45:55,1) = ones(11,1);
	v = 0.1 * v;
	
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
	
	result = @(t) [P * (C .* sin(sqrt(lambda) .* t + fi))];
	
	time  = 200;
	stime = 0;
	
	num_steps = time*10;
	time_step = time / num_steps;
   %{
	for t=0:num_steps-1
		data(1:100, t+1) = result(t*time_step + stime);
		
		h=plot(linspace(0, 100),result(t*time_step + stime));
		ylim([-1 1]);
        
		%saveas(h, strcat('plot', sprintf('%d', t), '.png'));
       pause(0.005);
        
    end
    %}
    %få noll vid vägg 2? data(1:101,num_steps+1)=zeros(101,1);
    
    %plotta 3d-plot med tiden som y axel vid diskreta tidpunkter
    %{
    hold on
    xlabel(['Partikel index']);
    ylabel(['$\frac{1}{\omega_o}$'], 'interpreter','latex');
    h_ylabel = get(gca,'YLabel');
set(h_ylabel,'FontSize',20); 
    zlabel(['Amplitud']);
    
    for t=0:time/100:time
        ett=ones(1,100);
        plot3(linspace(0, 100),t*ett,result(t*time_step + stime));
        zlim([-1 1]);
    end
    %}
    
    figure(2)
    
    for i = 1:100
        if sum(P(1, i)) < 0
            P(:, i) = -P(:, i);
        end
    end
    
    maxamp = P .* repmat(abs(C), 1, 100);
    hold on
    plot(repmat(1:100, 100, 1)', maxamp')
    
    
	%plot(repmat(linspace(stime, stime+time, num_steps), 100, 1)', data)
	
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