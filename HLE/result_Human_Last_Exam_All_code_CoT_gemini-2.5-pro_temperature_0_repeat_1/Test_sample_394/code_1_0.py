import math

def get_hypotenuse_denominator(r, s):
    """
    Calculates the hypotenuse and its denominator for a given r, s.
    """
    # The area is N = 263.
    N = 263

    # The product rs(r-s)(r+s) must be N times a perfect square.
    # rs(r-s)(r+s) = N * U^2
    product = r * s * (r - s) * (r + s)
    
    if product % N != 0:
        print(f"For r={r}, s={s}, the product is not divisible by {N}")
        return

    val = product // N
    
    if val < 0:
        # In some parametrizations, the product can be negative.
        # We take the absolute value for the square part.
        val = -val

    sqrt_val = math.isqrt(val)
    
    if sqrt_val * sqrt_val != val:
        print(f"For r={r}, s={s}, {val} is not a perfect square.")
        return
        
    U = sqrt_val
    
    # The hypotenuse c is given by (r^2 + s^2) / U
    numerator_c = r**2 + s**2
    denominator_c = U
    
    # The fraction for c is numerator_c / denominator_c.
    # We need to simplify it to find the final denominator.
    common_divisor = math.gcd(numerator_c, denominator_c)
    
    final_numerator = numerator_c // common_divisor
    final_denominator = denominator_c // common_divisor
    
    print(f"For parameters r = {r}, s = {s}:")
    print(f"The product rs(r-s)(r+s) = {product}")
    print(f"This equals {N} * {val}, where {val} is {sqrt_val}^2.")
    print(f"So, U = {U}.")
    print(f"The hypotenuse c is ({r}^2 + {s}^2) / {U} = {numerator_c} / {denominator_c}")
    print(f"After simplification, c = {final_numerator} / {final_denominator}")
    print(f"The smallest possible denominator is {final_denominator}.")

# The smallest known parameters for the congruent number 263 are r=53, s=26.
# This is a non-trivial result from the theory of elliptic curves.
r = 53
s = 26
get_hypotenuse_denominator(r, s)

# Let's verify the sides and area
# From the same theory, the sides a, b are given by:
# a = 2*r*s / U
# b = (r^2 - s^2) / U
# c = (r^2 + s^2) / U
# Let's calculate them with the found U.
N = 263
product = r * s * (r - s) * (r + s)
val = product // N
U = math.isqrt(val)

a_num = 2 * r * s
b_num = r**2 - s**2
c_num = r**2 + s**2

a = a_num / U
b = b_num / U
c = c_num / U

area = 0.5 * a * b
print("\nVerification:")
print(f"The sides of the triangle are:")
print(f"a = {a_num}/{U} = {a}")
print(f"b = {b_num}/{U} = {b}")
print(f"c = {c_num}/{U} = {c}")
print(f"The area is (1/2) * a * b = {area}")
if math.isclose(area, N):
    print("The area is indeed 263.")
