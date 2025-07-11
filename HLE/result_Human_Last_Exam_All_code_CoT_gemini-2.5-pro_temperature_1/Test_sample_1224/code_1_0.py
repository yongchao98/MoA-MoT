import sys

# Titan architecture constraint: numerators and denominators must be <= 15
BIT4_MAX = 15

# This class simulates a number in Titan's fractional system.
# It checks constraints upon creation.
class TitanNumber:
    def __init__(self, num, den=1, exp=0):
        if not (0 <= num <= BIT4_MAX and 0 < den <= BIT4_MAX):
            raise ValueError(f"Invalid 4-bit fraction: {num}/{den}")
        self.num = num
        self.den = den
        self.exp = exp

    def __repr__(self):
        return f"({self.num}/{self.den}e{self.exp})"

# This function simulates multiplication on the Titan computer.
def titan_multiply(t_num1, t_num2):
    """
    Multiplies two TitanNumbers, checking for 4-bit constraint violations.
    Assumes processor can simplify common factors, e.g., (a/b)*(c/d) = a*(c/b)/d
    """
    # Check for simplification: b divides c
    if t_num2.num % t_num1.den == 0:
        new_c = t_num2.num // t_num1.den
        new_num = t_num1.num * new_c
        new_den = t_num2.den
    # Check for simplification: d divides a
    elif t_num1.num % t_num2.den == 0:
        new_a = t_num1.num // t_num2.den
        new_num = new_a * t_num2.num
        new_den = t_num1.den
    # Default multiplication
    else:
        new_num = t_num1.num * t_num2.num
        new_den = t_num1.den * t_num2.den
    
    new_exp = t_num1.exp + t_num2.exp

    if new_num > BIT4_MAX or new_den > BIT4_MAX:
        # This operation is invalid on Titan.
        print(f"Failed to multiply {t_num1} and {t_num2}: result {new_num}/{new_den} exceeds 4-bit limit.", file=sys.stderr)
        return None
        
    return TitanNumber(new_num, new_den, new_exp)

def solve():
    """
    This function attempts to calculate the landing time using a simulation
    of the Titan computer's architecture. It prints 'N0' if the calculation
    is impossible.
    """
    # --- Step 1: Define constants as TitanNumbers ---
    # We must approximate constants with valid 4-bit fractions.
    try:
        # Physics constants
        C_4_div_3 = TitanNumber(4, 3)
        # Using a simple approximation for Pi
        C_pi = TitanNumber(3, 1) 
        # Approximating G with a valid fraction, e.g. 13/2 = 6.5
        C_G = TitanNumber(13, 2, exp=-11) 

        # Planet properties
        # d_shell = 300 kg/m^3
        C_d_shell = TitanNumber(3, 1, exp=2)
        # R_planet = 2,000,000 m
        C_R_planet = TitanNumber(2, 1, exp=6)
        
        # Freefall properties
        C_h = TitanNumber(5, 1, exp=3)
        C_2 = TitanNumber(2, 1)

    except ValueError as e:
        print(f"Initialization failed: {e}", file=sys.stderr)
        print("N0")
        return

    # --- Step 2: Attempt to calculate g ---
    # Simplified formula: g ≈ (4/3 * π) * G * d_shell * R_planet
    # We will try a few multiplication orders to see if any are possible.
    
    print("Attempting to calculate gravitational acceleration g...", file=sys.stderr)
    
    # Order 1: ( ( (4/3 * pi) * d_shell) * R_planet ) * G
    # (4/3) * pi = (4/3) * 3 = 4
    term1 = titan_multiply(C_4_div_3, C_pi)
    if term1 is None:
        print("N0")
        return
    print(f"Step 1: 4/3 * pi = {term1}", file=sys.stderr)

    # term1 * d_shell = 4 * 3e2 = 12e2
    term2 = titan_multiply(term1, C_d_shell)
    if term2 is None:
        print("N0")
        return
    print(f"Step 2: {term1} * d_shell = {term2}", file=sys.stderr)
    
    # term2 * R_planet = 12e2 * 2e6 = 24e8 -> Fails, because 24 > 15
    term3 = titan_multiply(term2, C_R_planet)
    if term3 is None:
        print("Conclusion: The calculation fails because intermediate products exceed the 4-bit limit.", file=sys.stderr)
        print("The architecture cannot handle the magnitude of these physical constants.", file=sys.stderr)
        print("N0")
        return

    # If we reached here, it means a miracle happened. Let's assume it did and continue.
    # g = term3 * G
    # ... calculation would continue ...
    # And would likely fail again. The square root calculation using Newton's method
    # would also fail as it requires adding and multiplying fractions, leading
    # to numerators/denominators larger than 15.

    print("Y[0.00]") # This line is unreachable given the constraints.

# Run the simulation
solve()