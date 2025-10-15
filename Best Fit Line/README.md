Here I tried to obtain a "best fit" line for the bivariate data, using a different approach.

Here "Best Fit" line is in the sense, that the line should be such that the sum of squares of perpendicular distances of the points from the line is minimized (rather than the classic regression approach).

My method is to, as if, physically rotate the data points about the origin and iteratively obtain the best fit vertical line that minimizes the $(x_i - x_{vert})^2$ for a given rotated data t. Now, keep on rotating the data points. Eventually, there will come a time when the "best fit" vertical line gives not only the $min(x_i - x_{vert})^2$ for the given rotated data t, but also the minimum most value for any rotated data set t.
