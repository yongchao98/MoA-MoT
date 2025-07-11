import math

# Define the maximum value for a 5-bit integer
MAX_VAL = 31

class TitanError(Exception):
    """Custom exception for Titan computer overflows."""
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

def titan_multiply(f1, f2):
    """
    Performs multiplication of two fractions (n, d) under Titan's 5-bit constraints.
    Raises TitanError on overflow.
    """
    n1, d1 = f1
    n2, d2 = f2

    # Check if input fractions are valid
    if not all(0 <= x <= MAX_VAL for x in [n1, d1, n2, d2]):
        raise TitanError(f"Input fraction contains value > {MAX_VAL}")

    # Simplify before multiplying to prevent intermediate overflow
    # g = gcd(a, b) -> a*c / b*d = (a/g)*(c) / (b/g)*d
    g1 = math.gcd(n1, d2)
    g2 = math.gcd(n2, d1)
    
    n1_s, d2_s = n1 // g1, d2 // g1
    n2_s, d1_s = n2 // g2, d1 // g2
    
    # Check if the simplified multiplication results in overflow
    num = n1_s * n2_s
    den = d1_s * d2_s

    if num > MAX_VAL or den > MAX_VAL:
        raise TitanError(
            f"Overflow multiplying ({n1}/{d1}) * ({n2}/{d2}). "
            f"Result {num}/{den} exceeds 5-bit limit."
        )
        
    return (num, den)

def titan_power(base, exp):
    """
    Performs exponentiation using Titan's multiplication rules.
    """
    if exp < 0:
        raise TitanError("Negative exponents are not supported.")
    if exp == 0:
        return (1, 1)
    
    result = base
    for i in range(exp - 1):
        # We print the operation that is being attempted
        print(f"Calculating step {i+2}: ({result[0]}/{result[1]}) * ({base[0]}/{base[1]})")
        result = titan_multiply(result, base)
    return result

def solve_pandora_problem():
    """
    Attempts to solve the Pandora gravity problem using Titan's rules.
    """
    # The mass calculation requires evaluating a term based on key ratios.
    # M_scaled = (rho_ratio - 1) + r_ratio^3
    rho_ratio_val = 4
    r_ratio_val = 20
    
    # We must represent these as Titan fractions
    rho_ratio = (rho_ratio_val, 1)
    r_ratio = (r_ratio_val, 1)
    
    print("The force calculation requires computing the planet's mass.")
    print("This in turn requires computing a term based on dimensionless ratios:")
    print(f"M_scaled = (rho_ratio - 1) + r_ratio^3")
    print("The numbers in this equation are:")
    print(f"rho_ratio = {rho_ratio[0]} / {rho_ratio[1]}")
    print(f"r_ratio = {r_ratio[0]} / {r_ratio[1]}")
    print("\nAttempting to calculate the term r_ratio^3 = (20/1)^3...")

    try:
        # This is the step that will fail. We try to compute 20^2 first.
        r_ratio_sq = titan_power(r_ratio, 2)
        # The code will not reach here, but if it did, it would continue.
        r_ratio_cubed = titan_multiply(r_ratio_sq, r_ratio)
        # ... and eventually calculate a force and an error.
        print("Calculation was successful (this should not happen).")

    except TitanError as e:
        print(f"\nCalculation failed: {e}")
        print("\nSince a key component of the mass formula cannot be computed due to hardware limitations, the gravitational force cannot be determined.")
        print("\n<<<N0>>>")

# Run the simulation
solve_pandora_problem()