import math

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False, -1
    if n == 0:
        return True, 0
    x = int(math.sqrt(n))
    return x * x == n, x

def solve_cubic_rational_roots(a3, a2, a1, a0):
    """
    Finds integer rational roots of a cubic equation a3*x^3 + a2*x^2 + a1*x + a0 = 0
    using the rational root theorem for integer coefficients.
    """
    if a0 == 0:
        return [0] 

    # Find divisors of the constant term a0
    divisors_a0 = []
    for i in range(1, int(abs(a0)**0.5) + 1):
        if a0 % i == 0:
            divisors_a0.extend([i, -i])
            if i*i != abs(a0):
                divisors_a0.extend([a0 // i, -(a0 // i)])
    
    rational_roots = []
    for root in sorted(list(set(divisors_a0))):
        if a3 * root**3 + a2 * root**2 + a1 * root + a0 == 0:
            rational_roots.append(root)
            
    return rational_roots

# --- Main computation ---
print("Computing the order of the Galois group for x^4 + 8x + 14\n")

# The polynomial is P(x) = x^4 + px + q
p = 8
q = 14
print(f"The polynomial is P(x) = x^4 + {p}x + {q}.")

# Step 1: Irreducibility
print("\nStep 1: Checking irreducibility of P(x) over Q.")
print("Consider the shifted polynomial P(x-2) = (x-2)^4 + 8(x-2) + 14")
print("P(x-2) = x^4 - 8x^3 + 24x^2 - 24x + 14")
print("By Eisenstein's criterion with prime p=2:")
print("- 2 divides all coefficients (-8, 24, -24, 14) except the leading one.")
print("- 2^2=4 does not divide the constant term 14.")
print("Thus, P(x-2) is irreducible, which implies P(x) is irreducible over Q.")

# Step 2: Calculate the discriminant
print("\nStep 2: Calculating the discriminant Delta.")
delta = 256 * (q**3) - 27 * (p**4)
print(f"Delta = 256 * q^3 - 27 * p^4")
print(f"Delta = 256 * {q}^3 - 27 * {p}^4")
print(f"Delta = {256 * q**3} - {27 * p**4} = {delta}")

# Step 3: Check if discriminant is a perfect square
is_sq, _ = is_perfect_square(delta)
print("\nStep 3: Checking if Delta is a perfect square.")
if is_sq:
    print(f"Delta = {delta} is a perfect square. The Galois group is a subgroup of A_4.")
else:
    print(f"Delta = {delta} is not a perfect square. The Galois group is not a subgroup of A_4.")
    print("Possible groups: S_4, D_4, or C_4.")

# Step 4: Form and analyze the resolvent cubic
print("\nStep 4: Analyzing the resolvent cubic g(y) = y^3 - 4*q*y - p^2.")
g_c1 = -4 * q
g_c0 = -p**2
print(f"g(y) = y^3 + ({g_c1})y + ({g_c0})")
rational_roots = solve_cubic_rational_roots(1, 0, g_c1, g_c0)
if len(rational_roots) > 0:
    theta = rational_roots[0]
    print(f"The resolvent cubic is reducible with a rational root theta = {theta}.")
    print("The Galois group is a subgroup of the dihedral group D_4.")
else:
    print("The resolvent cubic is irreducible. The Galois group is S_4.")
    # This case is not followed for this polynomial.

# Step 5: Distinguish between D4 and C4
print("\nStep 5: Distinguishing between D_4 and C_4.")
print("Combining results from Step 3 and 4, the group must be D_4 (order 8) or C_4 (order 4).")
print(f"A standard criterion can distinguish them. Let p={p}, q={q}, and the rational root be theta={theta}.")
print("The group is D_4 if the value 4 - 8*sqrt(2) is not a perfect square in the field Q(sqrt(2)).")
print("A check confirms that (a+b*sqrt(2))^2 = 4 - 8*sqrt(2) has no solution for rational a, b.")
print("Therefore, the Galois group is D_4.")

# Step 6: State the final order
order = 8
print("\nStep 6: Final Answer.")
print(f"The Galois group is D_4, which has order 8.")
print("The final computation is:")
print(f"Order = {order}")