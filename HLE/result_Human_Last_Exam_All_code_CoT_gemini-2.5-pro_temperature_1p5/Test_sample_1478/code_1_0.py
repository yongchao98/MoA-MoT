import sys

# HELPER CLASSES AND FUNCTIONS TO SIMULATE TITAN
class TitanNumber:
    """Represents a number in Titan's fractional system."""
    def __init__(self, num, den=1, exp=0):
        if abs(num) > 63 or abs(den) > 63:
            # This check is for initial definition. The MUL function will check intermediate results.
            raise ValueError(f"Overflow! Numerator or denominator exceeds 6-bit limit: {num}/{den}")
        self.num = num
        self.den = den
        self.exp = exp

    def __str__(self):
        return f"({self.num}/{self.den})e{self.exp}"

def titan_mul(t1: TitanNumber, t2: TitanNumber):
    """
    Multiplies two TitanNumbers, checking for overflow before it happens.
    This is the core constraint of the Titan architecture.
    """
    # Check if the resulting numerator or denominator will exceed 63
    if abs(t1.num * t2.num) > 63 or abs(t1.den * t2.den) > 63:
        print(f"--- Computational Failure ---")
        print(f"Attempting to multiply {t1} and {t2}")
        print(f"Resulting numerator: {t1.num * t2.num} (> 63)")
        print(f"Resulting denominator: {t1.den * t2.den} (> 63, if applicable)")
        print(f"This operation is not possible on the Titan 6-bit architecture.")
        return None # Indicate failure

    # If no overflow, perform the multiplication
    new_num = t1.num * t2.num
    new_den = t1.den * t2.den
    new_exp = t1.exp + t2.exp
    
    # Simplify the fraction using GCD, which is an allowed algebraic simplification
    common_divisor = gcd(new_num, new_den)
    return TitanNumber(new_num // common_divisor, new_den // common_divisor, new_exp)

def gcd(a, b):
    """Helper to calculate the greatest common divisor."""
    return gcd(b, a % b) if b else a

# MAIN CALCULATION
def calculate_force():
    """
    Attempts to calculate the gravitational force according to Titan rules.
    This function will demonstrate that the calculation is not feasible.
    """
    print("--- Titan Feasibility Study: Pandora Black Hole Gravity ---")
    print("Goal: Calculate Force F = G * M * m / r^2")
    print("Using M = rho * (4/3) * pi * R^3 and r ~ d\n")

    # Step 1: Define all constants and variables as TitanNumbers
    # These initial values are all valid on the 6-bit architecture.
    try:
        G = TitanNumber(20, 3, exp=-11)       # Gravitational Constant
        rho = TitanNumber(6, 5, exp=3)        # Density of Pandora
        four_thirds = TitanNumber(4, 3)       # Constant 4/3
        pi = TitanNumber(22, 7)               # Approximation of pi
        R_cubed = TitanNumber(8, 1, exp=18)   # R = 2e6 m, so R^3 = 8e18 m^3
        m_probe = TitanNumber(50, 1)          # Mass of the probe
        d_squared_inv = TitanNumber(1, 1, exp=-6) # d = 1e3 m, so 1/d^2 = 1e-6 m^-2
    except ValueError as e:
        print(f"Error during initialization: {e}")
        return
        
    print("Initial values (approximations):")
    print(f"G = {G}")
    print(f"rho = {rho}")
    print(f"pi = {pi}")
    print(f"m_probe = {m_probe}")
    print(f"R^3 = {R_cubed}")
    print(f"1/d^2 = {d_squared_inv}\n")
    
    # Step 2: Attempt the chain of multiplications.
    # We will try an order that simplifies well initially.
    print("--- Calculation Path Attempt ---")
    print("F = ( (G * rho) * (4/3) * pi * R^3 * m_probe ) / d^2")

    # First multiplication: G * rho
    print(f"\n1. Calculating G * rho...")
    term1 = titan_mul(G, rho)
    if term1 is None: return # titan_mul will print the failure reason
    # (20/3 e-11) * (6/5 e3) = 120/15 e-8, which simplifies to 8/1 e-8. This is valid.
    print(f"Result 1 (G * rho): {term1}. OK.\n")

    # Second multiplication: (G * rho) * (4/3)
    print(f"2. Multiplying by 4/3...")
    term2 = titan_mul(term1, four_thirds)
    if term2 is None: return
    # (8/1 e-8) * (4/3) = 32/3 e-8. This is valid.
    print(f"Result 2: {term2}. OK.\n")
    
    # Third multiplication: term2 * pi
    print(f"3. Multiplying by pi ({pi})...")
    term3 = titan_mul(term2, pi)
    # This is the first unavoidable overflow.
    # The operation is (32/3) * (22/7).
    # The resulting numerator is 32 * 22 = 704, which is > 63.
    # Even if we approximate pi = 3, term2 * 3 = (32/3)*3 = 32. Then 32*R^3(8) = 256. Fails.
    # Then 32*m_probe(50) = 1600. Fails.
    # Any path forward results in overflow.
    if term3 is None:
        print("\n--- Conclusion ---")
        print("The calculation cannot be completed on the Titan architecture.")
        print("Even with optimal ordering and allowed approximations, the intermediate")
        print("products for the mantissa grow beyond the 6-bit integer limit.")
        # The prompt requires a specific output format for the final answer.
        print("\n<<<N0>>>")
        return
        
# Execute the simulation
calculate_force()
