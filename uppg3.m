function uppg3a()
    figure(2)
    clf
	figure(1)
	clf
	hold on
	
	% particles for first calculation
	s = 20;
	% particles for last calculation
	N = 24;
    %number of eigen oscillations plotted
	num_plotted=6;
    
	for i=s:N
		[t, y] = calcStandingWaves(i);
        		
                
		% stretch all to the same length
		t = t * (N + 1)/(i + 1);
		
		% normalize maxumum amplitude for lowest frequency, so they will be easy to compare
		middle = y(1:(i+2), floor((i+2)/2));
		y = (y / max(middle)) * 0.9;
		
		% offset the waveforms, so they won't interfere
		y = y * 0.5 + i; 
		% 0.5 becuse we only want them going +-0.5, to allow proper y-axis numbering
        
        % ensure that we dont try to plot more eigen oscillation than we have
		to_plot=min(i,num_plotted);
        % plot to_plot first eigenoscillation
		plot(t(:,1:to_plot), y(:,1:to_plot));
	end
	
	hold off
	figure(2)
    hold on
	tau = linspace(0, 1, 1000);
    %number of particles gives number of eigen oscillations to plot
    num_particles=24;
    %first eigen oscillation we want to plot
    f=1
    %last eigen oscillation we want to plot
    l=num_plotted
    %plots eigen oscillations
	for i=f:l
    [t, y] = calcFrequencies(num_particles, i);
	plot(t,i+0.45*y);
    % +i to allow proper y-axis numbering
    % 0.45 to comfine plot to its number
   end

function [t, y] = calcStandingWaves(num_particles)
	A = triDiag(num_particles);
	
	[eigvec, ~] = eig(A);
	
	% create a list of arrays whose positions correspond to a node in the oscillations
	% 0 and num_particles+1 are endpoints
	t = 0:num_particles+1;
	% repeat for num_particles so we have them correspond to the number of eigenvectors
	% and transpose to make them ready for plot
	t = repmat(t, num_particles, 1)';
	
	% first node, 0
	y(1, 1:num_particles) = zeros(1, num_particles);
	% each eigenvector contains amplitudes for every node in a specific resonance oscillation,
	% add them so plot connects their amplitudes and connect them to the end nodes
	y(2:(num_particles+1), 1:num_particles) = eigvec;
	% last node, 0
	y((num_particles + 2), 1:num_particles) = zeros(1, num_particles);
	

function [t, y] = calcPositions(tau, num_particles, vector_index)
%plots the postions of individual particals over time for a specific
%eigen oscillation given by vector_index.
	assert(vector_index > 0, 'vector_index must be greater than zero');
	assert(vector_index <= num_particles, 'vector_index > num_particles, vector_index refers to a non-existant eigenvector');
	num_times = length(tau);
	
	A = triDiag(num_particles);
	
	%Hitta egenvärden:
	[eigvec, eigval] = eig(A);
	%matlab förutsätts ge normaliserade egenvektorerna
	
	%funktion som uttrycker positionen för massorna som avstånd m.a.p.
	%vänster
	%med tid uttryck num_particles 2pi/omega0
	%längdenhet l, l=L/(N+1) där N är antal
	%partiklar (X0(1) blir således 1, X0(2) 2, osv)
	
	% omega^2/omega_o^2 = eigval => omega/omega_o = sqrt(egival)
	eigfreqs = sqrt(eigval);
	F = eigfreqs * ones(num_particles, 1);
	%sin(tau)
	
	%tau = repmat(tau, num_particles, 1)
	
	% Only display one vector at a time
	eigvec = eigvec(1:num_particles, vector_index);
	F = F(vector_index, :);
	
	F   = repmat(F, 1, num_times);
	X0  = repmat((1:num_particles)', 1, num_times);
	
	X=@(t)(X0 + eigvec * (sin(2*pi*t .* F)));
	
	t = tau';
	y = X(tau)';
   
    
function [matrix] = triDiag(side_length)
	%Generera den tridiagonala matrisen:
	n = -ones(side_length - 1, 1);
	B = diag(n, 1);
	C = diag(n, -1);
	n = 2 * ones(side_length, 1);
	A = diag(n);
	matrix = A + B + C;
    

function [t,y] = calcFrequencies(num_particles,vector_index)
    K=triDiag(num_particles);
    [P D]=eig(K);
    freq=sqrt(D(vector_index,vector_index))*2*pi;
    %2*pi to plot in unit t=omega_o/2*pi
    t=linspace(0,10);
  
    y=sin(freq.*t);
        
    