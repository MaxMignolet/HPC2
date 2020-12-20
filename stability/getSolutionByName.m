% Get solution for the time step 'timestep'
% /!\ File must be located in directory 'results' and must be named
% 'u_X.dat' where 'X' is the timestep number.
%
function v = getSolutionByName(name, path)
	if nargin == 2
		fid = fopen([path name], 'r');
	else
		fid = fopen(['results/' name], 'r');
	end
    N = fread(fid,1,'int32');
    data = fread(fid,N*N,'double');
    fclose(fid);
    v = zeros(N,N);
    for a = 1:1:N
        for b = 1:1:N
            v(a,b) = data((a-1) * N + b);
        end
    end
    % if timeStep == 0
    %     v(find(v))
end
