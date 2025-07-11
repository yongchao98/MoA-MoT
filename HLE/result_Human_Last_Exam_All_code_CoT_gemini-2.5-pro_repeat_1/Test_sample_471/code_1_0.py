# The minimal number of critical points of a smooth function on a compact manifold
# is bounded below by the sum of its Betti numbers.
# For the 2-torus (T^2), the Betti numbers are well-known.

# b_0 is the number of connected components.
b0 = 1

# b_1 is the number of 1-dimensional or "circular" holes.
b1 = 2

# b_2 is the number of 2-dimensional "voids" or "cavities".
b2 = 1

# The sum of the Betti numbers gives the minimal number of critical points.
min_critical_points = b0 + b1 + b2

print("The minimal number of critical points is the sum of the Betti numbers: b0 + b1 + b2.")
print(f"For the 2-torus, the Betti numbers are b0 = {b0}, b1 = {b1}, and b2 = {b2}.")
print(f"The calculation is: {b0} + {b1} + {b2} = {min_critical_points}")