# A demonstration to show that monotonicity of f and g is not a sufficient
# condition for fp(f . g) = fp(f) ∩ fp(g).

# Define the poset L
# L is the set {0, 1, 2} with the usual ordering 0 <= 1 <= 2.
L = {0, 1, 2}

# Define the functions f and g.
# These functions are monotone on L.
# f is monotone: f(0)<=f(1) (0<=0), f(1)<=f(2) (0<=1)
# g is monotone: g(0)<=g(1) (1<=2), g(1)<=g(2) (2<=2)
def f(x):
    if x == 0: return 0
    if x == 1: return 0
    if x == 2: return 1

def g(x):
    if x == 0: return 1
    if x == 1: return 2
    if x == 2: return 2

# Define the composition (f . g)(x) = f(g(x))
def f_circ_g(x):
    return f(g(x))

# Helper function to find the set of fixed points for a given function
def find_fixed_points(func, domain):
    """Calculates the set of fixed points for a function."""
    return {x for x in domain if func(x) == x}

# 1. Calculate fp(f) and fp(g)
fp_f = find_fixed_points(f, L)
fp_g = find_fixed_points(g, L)

# 2. Calculate the intersection fp(f) ∩ fp(g)
intersection_fp = fp_f.intersection(fp_g)

# 3. Calculate fp(f . g)
fp_f_circ_g = find_fixed_points(f_circ_g, L)

# 4. Print the results and check for equality
print(f"Let L = {L} with the usual order.")
print(f"f(x) defined as: f(0)=0, f(1)=0, f(2)=1")
print(f"g(x) defined as: g(0)=1, g(1)=2, g(2)=2")
print("-" * 30)

# We want to check if fp(f . g) = fp(f) ∩ fp(g)
# The final equation is fp({0, 1}) = fp({0}) ∩ fp({2})
print(f"The fixed points of f are fp(f) = {fp_f}")
print(f"The fixed points of g are fp(g) = {fp_g}")
print(f"The intersection is fp(f) ∩ fp(g) = {intersection_fp}")
print(f"The fixed points of f.g are fp(f.g) = {fp_f_circ_g}")
print("-" * 30)
print("Final equation:")
print(f"fp(f.g) is the set of x such that f(g(x)) = x.")
print(f"f(g(0)) = f(1) = 0. So 0 is a fixed point.")
print(f"f(g(1)) = f(2) = 1. So 1 is a fixed point.")
print(f"f(g(2)) = f(2) = 1. So 2 is not a fixed point.")
print(f"fp(f.g) = {fp_f_circ_g}")

print(f"fp(f) is the set of x such that f(x)=x. So fp(f) = {fp_f}")
print(f"fp(g) is the set of x such that g(x)=x. So fp(g) = {fp_g}")

print(f"Resulting equation: {fp_f_circ_g} = {fp_f} ∩ {fp_g}")
print(f"Which simplifies to: {fp_f_circ_g} = {intersection_fp}")

if fp_f_circ_g == intersection_fp:
    print("\nThe equality holds.")
else:
    print("\nThe equality does NOT hold.")
    print("This demonstrates that 'f and g monotone' is not a sufficient condition.")