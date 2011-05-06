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
	
	xlabel('Antal partiklar', 'interpreter', 'latex');
	ylabel('$\omega$' ,'interpreter','latex');
	title(['egenfrekvenser f\"{o}r antal partiklar 1 till ', num2str(N)], 'interpreter','latex');
	
	
	% plots eigen oscillations
	for num_particles = 1:N
		% plot \omega/\omega_o for different numbers of particles
		for i = 1:num_particles
            % sets ticks 
			set(gca, 'XTick', 1:N)
            set(gca, 'YTick', 1:N)
			this_frequence = calcFrequencies(num_particles, i);
			plot(num_particles, this_frequence, '*');
		end
	end

	
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
	
function [frequence] = calcFrequencies(num_particles, vector_index)
    %returns sqrt of specified eigenvalue: \omega/\omega_o
	K     = triDiag(num_particles);
	[P D] = eig(K);
	% 2*pi to plot in unit t = omega_o / (2 * pi)
	frequence  = sqrt(D(vector_index, vector_index));