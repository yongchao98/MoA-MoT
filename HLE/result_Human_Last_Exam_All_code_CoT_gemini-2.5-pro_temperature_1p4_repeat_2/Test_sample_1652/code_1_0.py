import math

# Helper function to find the greatest common divisor for simplifying fractions
def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

# Helper functions for fraction arithmetic
def simplify_frac(n, d):
    """Simplifies a fraction (n/d)."""
    if d == 0:
        raise ZeroDivisionError
    common = gcd(n, d)
    return (n // common, d // common)

def add_frac(f1, f2):
    """Adds two fractions f1=(n1, d1) and f2=(n2, d2)."""
    n1, d1 = f1
    n2, d2 = f2
    n = n1 * d2 + n2 * d1
    d = d1 * d2
    return simplify_frac(n, d)

def mul_frac(f1, f2):
    """Multiplies two fractions."""
    n1, d1 = f1
    n2, d2 = f2
    n = n1 * n2
    d = d1 * d2
    return simplify_frac(n, d)

def div_frac(f1, f2):
    """Divides fraction f1 by f2."""
    n2, d2 = f2
    return mul_frac(f1, (d2, n2))

def isqrt(n):
    """
    Calculates the integer square root of a non-negative integer n.
    Uses an iterative method (Newton-Raphson) as sqrt is not available.
    """
    if n < 0:
        raise ValueError("isqrt() argument must be non-negative")
    if n == 0:
        return 0
    x = int(n**0.5) # Initial guess
    y = (x + n // x) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2
    return x

# Main logic to solve the physics problem
def solve():
    """
    Solves for the initial velocity u.
    The governing equation is Au^2 + Bu + C = 0
    """
    # Constants and approximations as fractions (n, d)
    # v = 5 m/s
    v_frac = (5, 1)
    # g = 9.8 m/s^2 = 98/10 = 49/5
    g_frac = (49, 5)
    # initial distance d = 300 m
    dist = 300
    # angle a = 60 degrees
    # sin(a) = sin(60) approx 13/15
    sin_a_frac = (13, 15)
    # sin(2a) = sin(120) approx 13/15
    sin_2a_frac = (13, 15)

    # Calculate coefficients A, B, C for Au^2 + Bu + C = 0
    # A = sin(2a)
    A_frac = sin_2a_frac
    # B = 2 * v * sin(a)
    B_frac = mul_frac((2, 1), mul_frac(v_frac, sin_a_frac))
    # C = -dist * g
    C_frac = mul_frac((-dist, 1), g_frac)

    print("The quadratic equation for the initial speed u is: Au^2 + Bu + C = 0")
    print(f"Calculated coefficients (as fractions n/d):")
    print(f"A = {A_frac[0]}/{A_frac[1]}")
    print(f"B = {B_frac[0]}/{B_frac[1]}")
    print(f"C = {C_frac[0]}/{C_frac[1]}")
    print("-" * 20)

    # Solve the quadratic equation u = (-B + sqrt(B^2 - 4AC)) / (2A)
    # Discriminant D = B^2 - 4AC
    b_sq_frac = mul_frac(B_frac, B_frac)
    four_ac_frac = mul_frac((4, 1), mul_frac(A_frac, C_frac))
    discriminant_frac = add_frac(b_sq_frac, (-four_ac_frac[0], four_ac_frac[1]))

    # Sqrt of discriminant
    d_n, d_d = discriminant_frac
    sqrt_d_n = isqrt(d_n)
    sqrt_d_d = isqrt(d_d)
    sqrt_discriminant_frac = (sqrt_d_n, sqrt_d_d)
    
    # Numerator: -B + sqrt(D)
    minus_b_frac = (-B_frac[0], B_frac[1])
    numerator_frac = add_frac(minus_b_frac, sqrt_discriminant_frac)

    # Denominator: 2A
    denominator_frac = mul_frac((2, 1), A_frac)

    # Final result for u
    u_frac = div_frac(numerator_frac, denominator_frac)
    u_val = u_frac[0] / u_frac[1]

    # Calculate memory usage z
    # 4 variables (A, B, C, u) of type frac. Each frac uses 6D.
    # struct frac { signed char (2D) n; unsigned char (2D) d; signed char (2D) e; } -> 6D
    z = 4 * 6
    
    # Round u to one decimal place
    u_rounded = round(u_val, 1)

    print(f"Result for u: {u_val:.4f} m/s")
    print(f"Rounded u: {u_rounded} m/s")
    print(f"Memory usage z: {z}D")
    
    return f"{u_rounded}:{z}"

# Run the solver and print the final answer
final_answer = solve()
print(f"\nFinal Answer (u:z) format:")
print(f"<<<{final_answer}>>>")
