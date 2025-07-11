# Plan:
# 1. Define a poset L. We'll use the set {0, 1, 2} with standard ordering.
# 2. Define two functions, f and g, which are both monotone on L.
# 3. Compute the fixed points for f, g, their intersection, and the fixed points of their composition f(g(x)).
# 4. Show that for these monotone functions, fp(f . g) is not equal to fp(f) ∩ fp(g).

# Let L be the set {0, 1, 2} with the usual order 0 < 1 < 2.
L = {0, 1, 2}

# Define a monotone function g. We choose g(0)=1, g(1)=2, g(2)=2.
# This is monotone because g(0)<=g(1) (1<=2) and g(1)<=g(2) (2<=2).
g_map = {0: 1, 1: 2, 2: 2}
def g(x):
    return g_map[x]

# Define a monotone function f. We choose f(0)=0, f(1)=0, f(2)=1.
# This is monotone because f(0)<=f(1) (0<=0) and f(1)<=f(2) (0<=1).
f_map = {0: 0, 1: 0, 2: 1}
def f(x):
    return f_map[x]

# Define the composition function h(x) = f(g(x))
def h(x):
    return f(g(x))

# Compute the set of fixed points for f
fp_f = {x for x in L if f(x) == x}

# Compute the set of fixed points for g
fp_g = {x for x in L if g(x) == x}

# Compute the intersection of the fixed point sets
intersection_fp = fp_f.intersection(fp_g)

# Compute the fixed points of the composition h = f . g
fp_h = {x for x in L if h(x) == x}

print("--- Counterexample for 'f and g monotone' ---")
print(f"Let poset L = {L} with the usual order.")
print(f"Let f be defined as: f(0)={f(0)}, f(1)={f(1)}, f(2)={f(2)}")
print(f"Let g be defined as: g(0)={g(0)}, g(1)={g(1)}, g(2)={g(2)}")
print("Both f and g are monotone functions.")
print("-" * 20)

print("The first set in the equation is fp(f) ∩ fp(g).")
print(f"Fixed points of f: fp(f) = {sorted(list(fp_f))}")
print(f"Fixed points of g: fp(g) = {sorted(list(fp_g))}")
# We need to output each number in the final sets for the equation
intersection_list = sorted(list(intersection_fp))
print(f"Equation Part 1: fp(f) ∩ fp(g) = {intersection_list}")
if not intersection_list:
    print("The set is empty.")
else:
    for item in intersection_list:
        print(f"  - {item}")
print("-" * 20)

print("The second set in the equation is fp(f . g).")
fp_h_list = sorted(list(fp_h))
print(f"Equation Part 2: fp(f . g) = {fp_h_list}")
if not fp_h_list:
    print("The set is empty.")
else:
    for item in fp_h_list:
        print(f"  - {item}")
print("-" * 20)

# Final comparison
is_equal = (fp_h == intersection_fp)
print(f"Result: Is fp(f . g) == fp(f) ∩ fp(g)? {is_equal}")
if not is_equal:
    print("\nThe equality does not hold. This demonstrates that 'f and g monotone' is not a sufficient condition.")