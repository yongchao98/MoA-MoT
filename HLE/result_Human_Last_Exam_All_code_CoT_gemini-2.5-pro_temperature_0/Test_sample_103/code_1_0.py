def find_fixed_points(func, domain):
    """Finds the fixed points of a function over a given domain."""
    return {x for x in domain if func(x) == x}

# Define the poset L = {0, 1}
L = {0, 1}

# Define two monotone functions f and g on L
# f(x) = 0 for all x
def f(x):
    return 0

# g(0) = 1, g(1) = 1
def g(x):
    if x == 0:
        return 1
    else:
        return 1

# Define the composition f . g
def f_dot_g(x):
    return f(g(x))

# Calculate the fixed point sets
fp_f = find_fixed_points(f, L)
fp_g = find_fixed_points(g, L)
fp_f_dot_g = find_fixed_points(f_dot_g, L)

# Calculate the intersection of fp(f) and fp(g)
intersection = fp_f.intersection(fp_g)

# Print the results and check the equality
print("This script tests the case where f and g are monotone.")
print(f"The poset is L = {L} with the usual order 0 <= 1.")
print("-" * 20)
print(f"f(x) = 0 for all x")
print(f"g(0) = 1, g(1) = 1")
print("-" * 20)
print(f"The set of fixed points of f is fp(f) = {fp_f}")
print(f"The set of fixed points of g is fp(g) = {fp_g}")
print(f"The intersection is fp(f) ∩ fp(g) = {intersection}")
print("-" * 20)
print(f"The composition (f.g)(x) has fixed points fp(f.g) = {fp_f_dot_g}")
print("-" * 20)
print("Checking the equality fp(f.g) = fp(f) ∩ fp(g):")
print(f"The equation becomes: {fp_f_dot_g} = {intersection}")
print(f"This is {fp_f_dot_g == intersection}.")
print("\nSince the equality is False, 'f and g monotone' is not a sufficient condition.")
