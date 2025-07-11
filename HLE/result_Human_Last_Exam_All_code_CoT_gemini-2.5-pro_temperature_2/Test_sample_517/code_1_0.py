import sympy

# Define the symbol n
n = sympy.symbols('n', real=True, positive=True)

# The probability p_n is asymptotically proportional to n^(-2/3).
# We can write p_n as a function f(n) for large n.
# f(n) = (2 / pi) * n**(-2/3)
# Here, we only need to compute the limit of this expression as n -> infinity.

# The expression for the approximate probability p_n
p_n_expression = (2 / sympy.pi) * n**(-sympy.Rational(2, 3))

# Calculate the limit as n tends to infinity
limit_p_n = sympy.limit(p_n_expression, n, sympy.oo)

# Print the step-by-step reasoning
print("Step 1: The random walk conditioned to avoid the origin has a uniform limiting angular distribution.")
print("Step 2: The path of the walk is approximated by a ray L_theta = (r*cos(theta), r*sin(theta)) where theta is uniform in [0, 2*pi).")
print("Step 3: The target is a disk of radius R_n = n^(1/3) centered at (n, 0).")
print("Step 4: The condition for the ray to intersect the disk is that the distance from (n,0) to the ray is less than R_n.")
print("The distance is d = n * |sin(theta)|.")
print("So we need n * |sin(theta)| <= n^(1/3), which means |sin(theta)| <= n^(-2/3).")
print("Step 5: For large n, this means the angle theta must be in a small interval of size proportional to n^(-2/3) around 0 or pi.")
print(f"Step 6: The probability p_n is the ratio of the favorable angle interval to 2*pi. So, p_n is proportional to n^(-2/3).")
print(f"p_n approx C * n^(-2/3) for some constant C. Here, C = 2/pi.")
print("Step 7: We calculate the limit of p_n as n -> infinity.")
print(f"lim_{{n->oo}} ({p_n_expression}) = {limit_p_n}")
