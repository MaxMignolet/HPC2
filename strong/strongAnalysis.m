dirPath = './outputs/';
filesNameArray = dir(dirPath);
nbFiles = length(filesNameArray);

% i = file index
% j = run index
% k = resolution scheme
j = ones(2, 1);

for i = 1:1:nbFiles
	fileName = filesNameArray(i).name;
	tmp = strsplit(fileName, '_');
	if length(tmp) == 3
		if tmp{2} == 'expl'
			k = 1;
		elseif tmp{2} == 'impl'
			k = 2;
		else
			k = -1; % should never happen in theory
		end
		file = fopen([dirPath fileName]);
		% 2nd line: nb processes
		% 3rd line: nb threads/process
		% 6th line: elapsed time
		% 7 lines in total

		while(true)
			line = fgetl(file); % 1
			if line == -1
				break;
			end
			line = fgetl(file); % 2
			nbProcess(k, j(k)) = str2double(line(22:end));
			line = fgetl(file); % 3
			nbThreads_per_process(k,  j(k)) = str2double(line(27:end));
			line = fgetl(file); % 4
			line = fgetl(file); % 5
			line = fgetl(file); % 6
			time(k, j(k)) = str2double(line(41:end-8));
			line = fgetl(file); % 7
			j(k) = j(k) + 1;
		end
	end
end

set(0,'defaultaxesfontsize', 15);
set(0,'defaulttextfontsize', 15);
set(0, 'defaultlegendfontsize', 15);
set(0,'defaultlinelinewidth', 2);

%% explicit part
k = 1; % resolution scheme

% graph of 1 process, 1:24 thread
% graph of 2 process, 1:24 thread
% graph of 4 process, 1:24 thread
% graph of 8 process, 1:24 thread

figure;
hold on
indices = find(nbProcess(k, :) == 1);
[sorted, I] = sort(nbThreads_per_process(k, indices));
plot(sorted, time(k, indices(I)), '*');
indices = find(nbProcess(k, :) == 2);
[sorted, I] = sort(nbThreads_per_process(k, indices));
plot(sorted, time(k, indices(I)), '*');
indices = find(nbProcess(k, :) == 4);
[sorted, I] = sort(nbThreads_per_process(k, indices));
plot(sorted, time(k, indices(I)), '*');
indices = find(nbProcess(k, :) == 8);
[sorted, I] = sort(nbThreads_per_process(k, indices));
plot(sorted, time(k, indices(I)), '*');

legend('1 process', '4 processes', '7 processes', '10 processes');
xlabel('Number of threads per process');
ylabel('Computation time [s]');
set(gca,'yscale','log')
grid
print(gcf,'scaling_expl1','-depsc')

% graph of 1:10 process, 3 thread
% graph of 1:10 process, 6 thread
% graph of 1:10 process, 12 thread
% graph of 1:10 process, 24 thread

figure;
hold on
indices = find(nbThreads_per_process(k, :) == 1);
[sorted, I] = sort(nbProcess(k, indices));
plot(sorted, time(k, indices(I)), '*');
indices = find(nbThreads_per_process(k, :) == 8);
[sorted, I] = sort(nbProcess(k, indices));
plot(sorted, time(k, indices(I)), '*');
indices = find(nbThreads_per_process(k, :) == 16);
[sorted, I] = sort(nbProcess(k, indices));
plot(sorted, time(k, indices(I)), '*');
indices = find(nbThreads_per_process(k, :) == 24);
[sorted, I] = sort(nbProcess(k, indices));
plot(sorted, time(k, indices(I)), '*');

legend('1 thread/process', '8 thread/process', '16 thread/process', '24 thread/process');
xlabel('Number of processes');
ylabel('Computation time [s]');
set(gca,'yscale','log')
grid
print(gcf,'scaling_expl2','-depsc')

%% implicit part
k = 2; % resolution scheme

% graph of 1 process, 1:24 thread
% graph of 2 process, 1:24 thread
% graph of 4 process, 1:24 thread
% graph of 8 process, 1:24 thread

figure;
hold on
indices = find(nbProcess(k, :) == 1);
[sorted, I] = sort(nbThreads_per_process(k, indices));
plot(sorted, time(k, indices(I)), '*');
indices = find(nbProcess(k, :) == 2);
[sorted, I] = sort(nbThreads_per_process(k, indices));
plot(sorted, time(k, indices(I)), '*');
indices = find(nbProcess(k, :) == 4);
[sorted, I] = sort(nbThreads_per_process(k, indices));
plot(sorted, time(k, indices(I)), '*');
indices = find(nbProcess(k, :) == 8);
[sorted, I] = sort(nbThreads_per_process(k, indices));
plot(sorted, time(k, indices(I)), '*');

legend('1 process', '4 processes', '7 processes', '10 processes');
xlabel('Number of threads per process');
ylabel('Computation time [s]');
set(gca,'yscale','log')
grid
print(gcf,'scaling_impl1','-depsc')

% graph of 1:10 process, 3 thread
% graph of 1:10 process, 6 thread
% graph of 1:10 process, 12 thread
% graph of 1:10 process, 24 thread

figure;
hold on
indices = find(nbThreads_per_process(k, :) == 1);
[sorted, I] = sort(nbProcess(k, indices));
plot(sorted, time(k, indices(I)), '*');
indices = find(nbThreads_per_process(k, :) == 8);
[sorted, I] = sort(nbProcess(k, indices));
plot(sorted, time(k, indices(I)), '*');
indices = find(nbThreads_per_process(k, :) == 16);
[sorted, I] = sort(nbProcess(k, indices));
plot(sorted, time(k, indices(I)), '*');
indices = find(nbThreads_per_process(k, :) == 24);
[sorted, I] = sort(nbProcess(k, indices));
plot(sorted, time(k, indices(I)), '*');

legend('1 thread/process', '8 thread/process', '16 thread/process', '24 thread/process');
xlabel('Number of processes');
ylabel('Computation time [s]');
set(gca,'yscale','log')
grid
print(gcf,'scaling_impl2','-depsc')
