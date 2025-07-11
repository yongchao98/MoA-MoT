import sympy as sp

# The problem asks for the minimal number of critical points for a smooth function f: T^2 -> R.
# We will solve this using a key result from Morse Theory and then verify it with a concrete example.

# --- Step 1: Use Morse Theory to find the lower bound ---
# Morse theory relates the number of critical points of a function on a manifold
# to the manifold's Betti numbers.
# The number of critical points is always greater than or equal to the sum of the Betti numbers.

# The Betti numbers for the 2-torus (T^2) are well-known:
# b_0: number of connected components
# b_1: number of 1-dimensional "handles" or loops
# b_2: number of 2-dimensional "voids"
b0 = 1
b1 = 2
b2 = 1

print("Step 1: Find the theoretical lower bound using Morse Theory.")
print("The Betti numbers for the 2-torus are:")
print(f"b_0 = {b0}")
print(f"b_1 = {b1}")
print(f"b_2 = {b2}")
print("-" * 40)

# The minimal number of critical points must be at least the sum of these numbers.
min_critical_points_bound = b0 + b1 + b2

print("The minimal number of critical points (N_crit) is bounded by the sum of Betti numbers:")
print("N_crit >= b_0 + b_1 + b_2")
print(f"N_crit >= {b0} + {b1} + {b2} = {min_critical_points_bound}")
print("-" * 40)


# --- Step 2: Show this minimum is achievable with an example ---
# To prove that 4 is the minimum, we must find a function that has exactly 4 critical points.
# Consider the function f(x, y) = cos(x) + cos(y) on the torus,
# represented by the square [0, 2*pi] x [0, 2*pi] with periodic boundaries.

print("Step 2: Construct a function that achieves this bound.")
print("Example function: f(x, y) = cos(x) + cos(y)")
print("We find its critical points in the domain [0, 2*pi) x [0, 2*pi).")
print("Critical points occur where the gradient is zero: [-sin(x), -sin(y)] = [0, 0]")

# Define symbolic variables and the function
x, y = sp.symbols('x y')
f = sp.cos(x) + sp.cos(y)

# The solutions for sin(x)=0 and sin(y)=0 in the domain are {0, pi}.
# This gives us 4 critical points.
critical_points = [(0, 0), (0, sp.pi), (sp.pi, 0), (sp.pi, sp.pi)]
print(f"\nThe critical points (x, y) are: {[(p[0], str(p[1])) for p in critical_points]}")

# Compute the Hessian matrix to classify the points
hessian_f = sp.hessian(f, (x, y))

# Classify each critical point by checking the eigenvalues of the Hessian
point_types = []
for p in critical_points:
    hess_at_p = hessian_f.subs({x: p[0], y: p[1]})
    eigenvals = list(hess_at_p.eigenvals().keys())
    if all(e > 0 for e in eigenvals):
        point_types.append("Local Minimum")
    elif all(e < 0 for e in eigenvals):
        point_types.append("Local Maximum")
    else:
        point_types.append("Saddle Point")

print("\nClassification of critical points:")
for p, ptype in zip(critical_points, point_types):
     print(f"  Point ({str(p[0]).ljust(3)}, {str(p[1]).ljust(3)}): {ptype}")

num_critical_points = len(critical_points)
print(f"\nThis function has a total of {num_critical_points} critical points.")
print("\nConclusion: Since theory proves the minimum is at least 4, and we found a function with exactly 4, the minimal number is 4.")
