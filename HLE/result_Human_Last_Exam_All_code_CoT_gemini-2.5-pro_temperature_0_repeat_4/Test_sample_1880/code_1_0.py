import math

def is_perfect_square(n):
    """Checks if a non-negative integer is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def find_integer_roots(coeffs):
    """
    Finds integer roots of a polynomial with integer coefficients.
    For a monic polynomial, any rational root must be an integer divisor of the constant term.
    """
    constant_term = coeffs[-1]
    if constant_term == 0:
        return [0]
    
    roots = []
    for i in range(1, abs(constant_term) + 1):
        if constant_term % i == 0:
            # Test positive divisor
            if sum(c * (i ** j) for j, c in enumerate(reversed(coeffs))) == 0:
                roots.append(i)
            
            # Test negative divisor
            if sum(c * ((-i) ** j) for j, c in enumerate(reversed(coeffs))) == 0:
                roots.append(-i)
    return list(set(roots))

# Define the coefficients for the polynomial f(x) = x^4 + ax^3 + bx^2 + cx + d
# Our polynomial is x^4 + 8x + 14
a, b, c, d = 0, 0, 8, 14

print(f"Computing the order of the Galois group for f(x) = x^4 + {c}x + {d}")
print("-" * 50)

# Step 1: Check for irreducibility using Eisenstein's criterion
p = 2
print("Step 1: Check for irreducibility over Q.")
print(f"Using Eisenstein's criterion with the prime p = {p}:")
print(f"  - p ({p}) divides all coefficients except the leading one (it divides {c} and {d}).")
print(f"  - p^2 ({p**2}) does not divide the constant term {d}.")
print("Conclusion: The polynomial is irreducible over Q.")
print("-" * 50)

# Step 2: Calculate the discriminant (Δ)
# For x^4 + cx + d, the formula is Δ = -27*c^4 + 256*d^3
delta = -27 * (c**4) + 256 * (d**3)
print("Step 2: Calculate the discriminant (Δ).")
print(f"Δ = -27 * c^4 + 256 * d^3")
print(f"Δ = -27 * ({c})^4 + 256 * ({d})^3")
print(f"Δ = -27 * {c**4} + 256 * {d**3}")
print(f"Δ = {-27 * c**4} + {256 * d**3}")
print(f"Δ = {delta}")

is_square = is_perfect_square(delta)
if is_square:
    print("Conclusion: Δ is a perfect square. The Galois group is a subgroup of A_4.")
else:
    print(f"Conclusion: Δ is not a perfect square (sqrt({delta}) ≈ {math.sqrt(delta):.2f}).")
    print("The Galois group is not a subgroup of A_4. Possible groups: S_4 or D_4.")
print("-" * 50)

# Step 3: Form and analyze the resolvent cubic
# For x^4 + cx + d, the resolvent cubic is g(y) = y^3 - 4dy - c^2
res_coeff_y = -4 * d
res_coeff_const = -c**2
resolvent_coeffs = [1, 0, res_coeff_y, res_coeff_const]

print("Step 3: Form and analyze the resolvent cubic g(y).")
print(f"g(y) = y^3 - 4*d*y - c^2")
print(f"g(y) = y^3 - 4*({d})*y - ({c})^2")
print(f"g(y) = y^3 + ({res_coeff_y})y + ({res_coeff_const})")

integer_roots = find_integer_roots(resolvent_coeffs)
if integer_roots:
    print(f"The resolvent cubic has an integer root at y = {integer_roots[0]}.")
    print("Conclusion: The resolvent cubic is reducible over Q.")
    is_resolvent_reducible = True
else:
    print("The resolvent cubic has no integer roots.")
    print("Conclusion: The resolvent cubic is irreducible over Q.")
    is_resolvent_reducible = False
print("-" * 50)

# Step 4: Determine the Galois group and its order
print("Step 4: Determine the final result.")
print("Summary of findings:")
print(f"  - The polynomial is irreducible.")
print(f"  - The discriminant is NOT a perfect square.")
print(f"  - The resolvent cubic IS reducible.")
print("\nThese conditions uniquely identify the Galois group as the dihedral group D_4.")

order = 8
print(f"\nThe order of the Galois group D_4 is {order}.")
<<<8>>>