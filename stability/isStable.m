dirPath = './results3/';
filesNameArray = dir(dirPath);

nbFiles = length(filesNameArray);
% is_stable = zeros(nbFiles-2, 1);
% N = zeros(nbFiles-2, 1, 'uint32');
% CourantNumber = zeros(nbFiles-2, 1, 'double');

j = 1;
for i = 1:1:nbFiles
	fileName = filesNameArray(i).name;
	tmp = strsplit(fileName, '_');
	if length(tmp) > 2
		N(j) = str2double(tmp{2});
		CourantNumber(j) = str2double(tmp{3});
		sol = getSolutionByName(fileName, dirPath);
		if any(isnan(sol)) == 1
			is_stable(j) = 0;
		else
			is_stable(j) = 1;
		end
		j = j + 1;
	end
end

% display(is_stable)
stable_indices = find(is_stable);
unstable_indices = find(~is_stable);

stable_points = zeros(2, length(stable_indices));
unstable_points = zeros(2, length(unstable_indices));

stable_points(1, :) = N(stable_indices);
stable_points(2, :) = CourantNumber(stable_indices);
unstable_points(1, :) = N(unstable_indices);
unstable_points(2, :) = CourantNumber(unstable_indices);

set(0,'defaultaxesfontsize', 15);
set(0,'defaulttextfontsize', 15);
set(0,'defaultlinelinewidth', 2);

figure;
scatter(stable_points(1, :), stable_points(2, :), [], 'g', 'filled');
hold on
scatter(unstable_points(1, :), unstable_points(2, :), [], 'r', 'filled');
set(gca,'yscale','log')
grid

xlabel('Number of elements')
ylabel('Courant number')
% axis([0 220. 5e-3 2])
print(gcf,'stability3','-depsc')
