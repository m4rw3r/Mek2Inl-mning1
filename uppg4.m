function apa()
	clf
	figure(1)
	clf
	
	% Number of particles
	N = 100;
	
	% Generera tridiagonala matrisen
	A = triDiag(N);
	
	% Homogenledet
	noll = zeros(N,1);
	% Hastighetsledet
	v    = zeros(N,1);
	% Ger mitterta partiklarna en hastighet vid t = 0
	v(45:55,1) = ones(11,1);
	v          = 0.1 * v;
	
	%__________
	% Vi har ekvationerna:
	% P C .* sin(fi) = noll                        (1)
	% P sqrt(lambda) * omega0 .* C .* cos(Fi) = v  (2)
	% (1) => sin(fi) = 0 (ty C != 0) => fi = pi*n
	% (2) ger:
	% P * C .* cos(Fi) = v ./ (sqrt(lambda) * omega0)
	% C .* cos(Fi) = P \ (v ./ (sqrt(lambda) * omega0))
	% Varpå C fås från insättning av Fi, från (1), och Fi fås från C.
	%___________
	
	rhs(1:N, 1) = noll;
	rhs(1:N, 2) = v;
	
	[P D]  = eig(A);
	lambda = D * ones(N, 1);
	
	solution = P \ rhs;
	
	% Vi låter omega0 vara tidsenhet, (fås hastighet i enhet meter*omega0)
	k = sqrt(lambda); % * omega0
	C = solution(:, 2) ./ k;
	fi = 0;
	
	result = @(t) [P * (C .* sin(sqrt(lambda) .* t + fi))];
	
	time      = 200;
	stime     = 0;
	num_steps = time * 10;
	time_step = time / num_steps;
	
	for t = 0:(num_steps - 1)
		data(1:N, t + 1) = result(t * time_step + stime);
		
		h = plot((1:N)', result(t * time_step + stime));
		ylim([-1 1]);
		
		%saveas(h, strcat('plot', sprintf('%d', t), '.png'));
		pause(0.005);
	end
	%}
	% Få noll vid vägg 2? data(1:101, num_steps + 1) = zeros(101, 1);
	
	% Plotta 3d-plot med tiden som y axel vid diskreta tidpunkter
	%{
	clf
	hold on
	xlabel(['Partikel-index']);
	ylabel(['$\frac{1}{\omega_o}$'], 'interpreter', 'latex');
	h_ylabel = get(gca, 'YLabel');
	set(h_ylabel, 'FontSize', 20); 
	zlabel(['Amplitud']);
	
	for t = 0:time/N:time
		ett = ones(1, N);
		plot3(linspace(0, N), t * ett, result(t * time_step + stime));
		zlim([-1 1]);
	end
	%}
	
	figure(2)
	
	% Force the first eigen-oscillation to be positive,
	% matlab might not give them all the same sign
	for i = 1:N
		if sum(P(1, i)) > 0
			P(:, i) = -P(:, i);
		end
	end
	
	maxamp = P .* repmat(C, 1, N);
	hold on
	plot(repmat(1:N, N, 1)', maxamp')
	
function [matrix] = triDiag(side_length)
	% Generera den tridiagonala matrisen:
	n = -ones(side_length - 1, 1);
	B = diag(n, 1);
	C = diag(n, -1);
	n = 2 * ones(side_length, 1);
	A = diag(n);
	matrix = A + B + C;