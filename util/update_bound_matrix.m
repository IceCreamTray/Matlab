function update_bound_matrix(Min, xdim, ydim, x, y, value)
	Min(x, y) = value;
	update_bound_matrix_bounds(Min, xdim, ydim);
end