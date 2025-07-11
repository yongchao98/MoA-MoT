import math

def is_perfect_square(n):
    """Checks if a non-negative integer is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

# The polynomial is f(x) = x^4 + 8x + 14.
# In the form x^4 + ax^3 + bx^2 + cx + d, the coefficients are:
a = 0
b = 0
c = 8
d = 14

print(f"Analyzing the polynomial f(x) = x^4 + {c}x + {d}\n")

# Step 1: Check for irreducibility
print("--- Step 1: Check for Irreducibility ---")
# We apply Eisenstein's criterion with the prime p=2.
# 1. p=2 does not divide the leading coefficient (1).
# 2. p=2 divides all other coefficients (c=8, d=14).
# 3. p^2=4 does not divide the constant term (d=14).
print("By Eisenstein's criterion with p=2, the polynomial is irreducible over Q.")
print("This means the Galois group is a transitive subgroup of S_4.\n")


# Step 2: Compute the discriminant
print("--- Step 2: Compute the Discriminant ---")
# For a depressed quartic x^4 + cx + d, the discriminant is Δ = 256*d^3 - 27*c^4.
discriminant = 256 * (d**3) - 27 * (c**4)
print(f"The discriminant formula is Δ = 256*d^3 - 27*c^4")
print(f"Δ = 256 * ({d}^3) - 27 * ({c}^4)")
print(f"Δ = 256 * {d**3} - 27 * {c**4}")
print(f"Δ = {256 * d**3} - {27 * c**4}")
print(f"Δ = {discriminant}")

is_square = is_perfect_square(discriminant)
print(f"Is the discriminant a perfect square in Q? {is_square}")
if not is_square:
    print("Since Δ is not a perfect square, the Galois group is not a subgroup of A_4.")
    print("Possible groups: S_4, D_4, C_4.\n")
else:
    print("Since Δ is a perfect square, the Galois group is a subgroup of A_4.")
    print("Possible groups: A_4, V_4.\n")


# Step 3: Form the resolvent cubic
print("--- Step 3: Form the Resolvent Cubic ---")
# For x^4 + cx + d, the resolvent cubic is g(y) = y^3 - 4*d*y - c^2.
res_c0 = -(c**2)
res_c1 = -4 * d
print(f"The resolvent cubic formula is g(y) = y^3 - 4*d*y - c^2")
print(f"g(y) = y^3 - 4*({d})*y - ({c}^2)")
print(f"g(y) = y^3 + {res_c1}y + {res_c0}\n")


# Step 4: Analyze the resolvent cubic
print("--- Step 4: Analyze the Resolvent Cubic ---")
# Check for rational roots using the Rational Root Theorem.
# Possible rational roots must be integer divisors of the constant term -64.
constant_term = abs(res_c0)
divisors = [i for i in range(1, constant_term + 1) if constant_term % i == 0]
rational_roots = []
for r in divisors + [-i for i in divisors]:
    if r**3 + res_c1 * r + res_c0 == 0:
        rational_roots.append(r)

if len(rational_roots) > 0:
    is_resolvent_reducible = True
    print(f"The resolvent cubic has rational root(s): {rational_roots}.")
    print("Therefore, the resolvent cubic is reducible over Q.\n")
else:
    is_resolvent_reducible = False
    print("The resolvent cubic has no rational roots, so it is irreducible over Q.\n")


# Step 5: Determine the Galois group and its order
print("--- Step 5: Determine the Galois Group and its Order ---")
print("Summary of findings:")
print(f"  - The polynomial is irreducible.")
print(f"  - The discriminant ({discriminant}) is NOT a perfect square.")
print(f"  - The resolvent cubic is REDUCIBLE over Q.")

# Classification:
# - Irreducible, Δ not square, resolvent irreducible -> S_4 (order 24)
# - Irreducible, Δ not square, resolvent reducible   -> D_4 or C_4 (order 8 or 4)
# - Irreducible, Δ square, resolvent irreducible     -> A_4 (order 12)
# - Irreducible, Δ square, resolvent reducible       -> V_4 (order 4)
#
# To distinguish D_4 from C_4, we check if all roots of the resolvent are rational.
# If only one root is rational, the group is D_4.
# Our resolvent is g(y) = y^3 - 56y - 64. We found the rational root y=8.
# Factoring gives (y-8)(y^2 + 8y + 8).
# The roots of y^2 + 8y + 8 are (-8 ± sqrt(64-32))/2 = -4 ± 2*sqrt(2), which are not rational.
# Thus, there is only one rational root.
galois_group = "D_4"
order = 8

print(f"\nThis combination of properties uniquely identifies the Galois group as the dihedral group {galois_group}.")
print(f"The order of the Galois group {galois_group} is {order}.")
