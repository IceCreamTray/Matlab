function update_bound_matrix(Min, xdim, ydim, x, y, value)
	Min(x, y) = value;
	
	doUpdate = 0;
	nx = x;
	ny = y;
	
	if (x <= 2)
		doUpdate = 1;
		
		nx = (x + (xdim - 2));
	elseif (x >= (xdim - 1))
		doUpdate = 1;
		
		nx = (2 - (xdim - x));
	end
	
	if (y <= 2)
		doUpdate = 1;
		
		ny = (y + (ydim - 2));
	elseif (y >= (ydim - 1))
		doUpdate = 1;
		
		ny = (2 - (ydim - y));
	end
	
	if (doUpdate == 1)
		Min(nx, ny) = value;
	end
	
	%update_bound_matrix_bounds(Min, xdim, ydim);
end