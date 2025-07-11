import math

# A helper function to represent a fraction and check constraints.
def create_fraction(numerator, denominator, name=""):
    """Creates a fraction, checking if it respects Titan's 5-bit constraints."""
    if not (0 <= numerator <= 31 and 0 < denominator <= 31):
        raise ValueError(f"Fraction {name} ({numerator}/{denominator}) violates 5-bit constraint.")
    return (numerator, denominator)

# A helper for the greatest common divisor, used for simplification.
def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

# A helper function for Titan multiplication.
def titan_multiply(frac1, frac2):
    """
    Multiplies two fractions, simplifying before multiplication to avoid overflow.
    This simulates how the Titan computer would handle the operation.
    """
    n1, d1 = frac1
    n2, d2 = frac2

    # Simplify before multiplying to stay within 5-bit constraints.
    # Simplify n1 with d2
    common1 = gcd(n1, d2)
    n1_s = n1 // common1
    d2_s = d2 // common1

    # Simplify n2 with d1
    common2 = gcd(n2, d1)
    n2_s = n2 // common2
    d1_s = d1 // common2
    
    # Perform the multiplication on the simplified terms
    res_n = n1_s * n2_s
    res_d = d1_s * d2_s
    
    # Check if the result is valid
    create_fraction(res_n, res_d)
    
    return (res_n, res_d)

def solve():
    """
    Solves the problem using the Titan computational model.
    """
    print("Yes, the Titan computer can solve this problem.")
    print("--- Physics Formulation ---")
    
    # Derivation
    # Time to fall h: t = sqrt(2*h/g)
    # Distance travelled: d = (1/2)*a*t^2 = (1/2)*(F/m)*t^2
    # Substitute t^2: d = (1/2)*(F/m)*(2*h/g) = F*h/(m*g)
    # Solve for F: F = d*m*g/h
    # Mass m = density * volume = rho * (4/3)*pi*r^3
    # With d=20, h=10 -> d/h = 2. So F = 2*m*g
    # With rho=0.9kg/cm^3 = 9e5 kg/m^3, r=0.5cm=1/200m -> m = 9e5 * (4/3)*pi*(1/200)^3 = (3/20)*pi
    # So, F = 2 * (3/20)*pi * g = (3/10)*pi*g

    print("The governing equation for the force is: F = (3/10) * pi * g\n")

    print("--- Titan Calculation ---")
    print("Step 1: Represent all constants as 5-bit integer fractions.")
    
    f_3_10 = create_fraction(3, 10, "3/10")
    
    # We must approximate pi and g. To ensure the calculation is possible,
    # we select values that allow for simplification.
    # pi approx. 28/9 = 3.111...
    # g approx. 10/1 = 10.0
    pi_approx = create_fraction(28, 9, "pi")
    g_approx = create_fraction(10, 1, "g")

    print(f"Using approximation pi ≈ {pi_approx[0]}/{pi_approx[1]}")
    print(f"Using approximation g ≈ {g_approx[0]}/{g_approx[1]}\n")

    # The calculation is F = (3/10) * (28/9) * (10/1)
    print("Step 2: Perform multiplication, simplifying at each stage to respect constraints.")
    print("Calculation: (3/10) * (10/1) * (28/9)")
    
    # First, multiply (3/10) by (10/1) for a convenient simplification
    term1 = f_3_10
    term2 = g_approx
    intermediate_result = titan_multiply(term1, term2)
    print(f"  ({term1[0]}/{term1[1]}) * ({term2[0]}/{term2[1]}) simplifies to {intermediate_result[0]}/{intermediate_result[1]}")

    # Next, multiply the result by our pi approximation
    term3 = pi_approx
    final_result = titan_multiply(intermediate_result, term3)
    print(f"  ({intermediate_result[0]}/{intermediate_result[1]}) * ({term3[0]}/{term3[1]}) simplifies to {final_result[0]}/{final_result[1]}\n")
    
    print("--- Final Result ---")
    print("The final equation with the chosen fractional approximations is:")
    final_equation_str = (
        f"{f_3_10[0]} / {f_3_10[1]} * "
        f"{pi_approx[0]} / {pi_approx[1]} * "
        f"{g_approx[0]} / {g_approx[1]} = "
        f"{final_result[0]} / {final_result[1]}"
    )
    print(final_equation_str)

    # For the final answer, calculate the error
    f_titan = final_result[0] / final_result[1]
    f_true = (3/10) * math.pi * 9.8
    abs_error = abs(f_titan - f_true)
    
    # This part is for the final answer format and not part of the simulation
    # print(f"\nSmallest absolute error e = {abs_error}")

solve()