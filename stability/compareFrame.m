function residuals = compareFrame(x, y)
	N = size(x, 1);
	diff_xy = 0;
	norm_x = 0;
	norm_y = 0;
	for i = 1:N
		for j = 1:N
			diff_xy = diff_xy + (x(i, j) - y(i, j))^2;
			norm_x = norm_x + x(i, j)^2;
			norm_y = norm_y + y(i, j)^2;
		end
	end
	residuals = diff_xy / ((norm_x + norm_y) / 2);
end
