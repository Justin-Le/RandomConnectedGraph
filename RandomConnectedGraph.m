%% Example

clc; clearvars; close all;

seed = 0;
n = 20;
sparsity = 0.4; % must be a positive real number strictly less than 0.5

[gr, Laplacian] = generateRandomConnectedGraph(seed, n, sparsity);

plot(gr)
eigvals = eig(Laplacian);
fprintf('Connectivity (smallest positive eigenvalue of the Laplacian): %d.\n\n', eigvals(2));

function [gr, Laplacian] = generateRandomConnectedGraph(seed, n, sparsity)
    % Generates a random undirected connected graph with a user-specified
    % sparsity level.
    
    % Inputs: the RNG seed; the number of nodes; 
    % the sparsity level (a positive real number strictly less than 0.5)
    
    % Outputs: the Matlab graph object; the graph Laplacian matrix.
    
	rng(seed);

	adj = round(rand(n) - sparsity);
	adj = triu(adj) + triu(adj, 1)';
	adj = adj - diag(diag(adj));
	gr = graph(adj);
	bins = conncomp(gr);
	num_components = max(bins);

	if num_components > 1
	    new_neighbors = zeros(1, num_components);

	    for k = 1:num_components
	        component_indices = find(bins == k);

	        if isscalar(component_indices)
	            new_neighbors(k) = component_indices;
	        else
	            new_neighbors(k) = randsample(component_indices, 1);
	        end
	    end
	end

	for k = 1:num_components-1
	    adj(new_neighbors(k), new_neighbors(k+1)) = 1;
	    adj(new_neighbors(k+1), new_neighbors(k)) = 1;
	end

	gr = graph(adj);
	Laplacian = diag(sum(adj, 2)) - adj;
end
