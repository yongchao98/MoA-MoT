# The Betti numbers for the 2-torus (T^2)
# b0: Number of connected components
b0 = 1
# b1: Number of 1-dimensional holes (tunnels)
b1 = 2
# b2: Number of 2-dimensional holes (voids)
b2 = 1

# According to the Morse inequalities, the number of critical points of each index
# (c0 for minima, c1 for saddles, c2 for maxima) must be at least the
# corresponding Betti number.
# c0 >= b0
# c1 >= b1
# c2 >= b2
# Therefore, the minimal total number of critical points is the sum of the Betti numbers.
min_minima = b0
min_saddles = b1
min_maxima = b2

total_critical_points = min_minima + min_saddles + min_maxima

# The final equation shows the minimal number of each type of critical point
# and the total minimal number of critical points.
print(f"The minimal number of critical points for a smooth function on a 2-torus is given by the sum of its Betti numbers.")
print(f"Minimal number of minima (index 0) = {min_minima}")
print(f"Minimal number of saddles (index 1) = {min_saddles}")
print(f"Minimal number of maxima (index 2) = {min_maxima}")
print(f"Total minimal number of critical points:")
print(f"{min_minima} + {min_saddles} + {min_maxima} = {total_critical_points}")