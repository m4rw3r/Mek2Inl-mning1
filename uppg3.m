function uppg3a()
	figure(2)
	clf
	figure(1)
	clf
	
	% particles for first calculation
	s = 20;
	% particles for last calculation
	N = 24;
	%number of eigen oscillations plotted
	num_plotted = 6;
	
	%set titles and labels for figure 1
	title([num2str(num_plotted), ' f\"{o}rsta egensv\"{a}ningarna som st\aa{}ende v\aa{}gor f\"{o}r ', num2str(s), ' till ', num2str(N), ' partiklar'], 'interpreter', 'latex');
	xlabel('Avst\aa{}nd fr\aa{}n v\"{a}nsterv\"{a}ggen', 'interpreter', 'latex');
	ylabel('Olika antal partiklar samt deras amplituder');
	
	for i = s:N
		[t, y] = calcStandingWaves(i);
		
		% normalize maxumum amplitude for lowest frequency, so they will be easy to compare
		middle = y(1:(i + 2), floor((i + 2) / 2));
		y = (y / max(middle)) * 0.9;
		
		% ensure that we dont try to plot more eigen oscillation than we have
		to_plot = min(i,num_plotted);
		
		subplot(N - s + 1, 1 , i + 1 - s)
		hold on
		
		ylabel([num2str(i), ' par.']);
		set(gca,'XTick',1:N)
		plot(t(:, 1:to_plot), y(:, 1:to_plot), '-*');
	end
	
	xlabel('Partikel-index')
	
	hold off
	figure(2)
	
	hold on;
	
	first = 1;
	last  = num_plotted;
	
	xlabel('tid, $\frac{2\pi}{\omega_o}$', 'interpreter', 'latex');
	ylabel('index f\"{o}r egenfrekvens samt normaliserad amplitud' ,'interpreter','latex');
	title(['egenfrekvenser ', num2str(first), ' till ', num2str(last) ,' f\"{o}r antal partiklar ', num2str(s), ' till ', num2str(N), ' i olika f\"{a}rger'], 'interpreter','latex');
	
	% to plot differenet frequencies with different colors.
	colors = ['b', 'r', 'g', 'k', 'c', 'y', 'm'];
	
	% plots eigen oscillations
	for num_particles = s:N
		% to plot different numbers of particles in different subplots
		subplot(N - s + 1, 1, num_particles + 1 - s)
		ylabel([num2str(num_particles), ' par.']);
		hold on
		
		for i = first:last
			colorindex = mod(i - last, length(colors)) + 1;
			
			% sets ticks 
			set(gca, 'XTick', 0:10)
			[t, y] = calcFrequencies(num_particles, i);
			plot(t, y, colors(colorindex));
		end
	end
	
	% sets common xlabel
	 xlabel('$\frac{2\pi}{\omega_o}$', 'interpreter','latex');
	
function [t, y] = calcStandingWaves(num_particles)
	A = triDiag(num_particles);
	
	[eigvec, ~] = eig(A);
	
	t = 0:num_particles + 1;
	t = repmat(t, num_particles, 1)';
	
	y(1, 1:num_particles) = zeros(1, num_particles);
	% each eigenvector contains amplitudes for every node in a specific resonance oscillation,
	% add them so plot connects their amplitudes and connect them to the end nodes
	y(2:(num_particles+1), 1:num_particles) = eigvec;
	y((num_particles + 2), 1:num_particles) = zeros(1, num_particles);
	
function [m] = triDiag(side_length)
	%Generera den tridiagonala matrisen:
	n = -ones(side_length - 1, 1);
	B = diag(n, 1);
	C = diag(n, -1);
	n = 2 * ones(side_length, 1);
	A = diag(n);
	m = A + B + C;
	
function [t, y] = calcFrequencies(num_particles, vector_index)
	K     = triDiag(num_particles);
	[P D] = eig(K);
	% 2*pi to plot in unit t = omega_o / (2 * pi)
	freq  = sqrt(D(vector_index, vector_index)) * 2 * pi;
	t     = linspace(0, 10);
	y     = sin(freq .* t);