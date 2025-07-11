import math

class TitanError(Exception):
    """Custom exception for Titan architecture violations."""
    pass

class TitanFrac:
    """
    Represents a number in Titan's fractional system: num/den * 10^exp.
    Enforces the 4-bit constraint on numerators and denominators.
    """
    def __init__(self, num, den, exp=0):
        if not (0 <= num <= 15 and 1 <= den <= 15):
            raise TitanError(f"Initialization failed: Numerator ({num}) or denominator ({den}) out of 4-bit range [0, 15].")
        self.num = num
        self.den = den
        self.exp = exp

    def __mul__(self, other):
        """
        Multiplies two TitanFrac numbers and checks for 4-bit constraint violation.
        """
        new_num = self.num * other.num
        new_den = self.den * other.den
        
        # In a real system, one would simplify using Greatest Common Divisor.
        common_divisor = math.gcd(new_num, new_den)
        simplified_num = new_num // common_divisor
        simplified_den = new_den // common_divisor

        if simplified_num > 15 or simplified_den > 15:
            # According to the rules, this operation is invalid. The architecture cannot proceed.
            raise TitanError(f"FATAL: Product ({self.num}/{self.den})*({other.num}/{other.den}) -> {simplified_num}/{simplified_den} exceeds 4-bit limit.")
        
        new_exp = self.exp + other.exp
        return TitanFrac(simplified_num, simplified_den, new_exp)

    def __repr__(self):
        """Formats the number for printing in an equation."""
        if self.exp == 0:
            return f"({self.num}/{self.den})"
        else:
            return f"({self.num}/{self.den} * 10^{self.exp})"

def solve():
    """
    Attempts to calculate Pandora's escape velocity using the Titan architecture simulation.
    """
    print("This script simulates the Titan 4-bit architecture to determine if the escape velocity of Pandora can be calculated.")
    print("The escape velocity squared is given by v_e^2 = 2*G*M/R.")
    print("We simplify the calculation by assuming the planet's mass M is dominated by its shell (M_core is negligible).")
    print("The formula becomes: v_e^2 ≈ (8/3) * G * pi * rho_shell * R^2")
    print("-" * 50)
    
    try:
        # Define constants using the most viable fractional approximations that fit in 4 bits
        # G ≈ 6.674e-11 m^3/kg/s^2. Best 4-bit approximation is 13/2 = 6.5
        G = TitanFrac(13, 2, -11)
        # pi ≈ 3.14159. The simplest integer approximation is 3.
        PI = TitanFrac(3, 1)
        # rho_shell = 300 kg/m^3
        RHO_SHELL = TitanFrac(3, 1, 2)
        # R = 2e6 m, so R^2 = 4e12 m^2
        R_SQUARED = TitanFrac(4, 1, 12)
        # The constant 8/3 from the volume formula (4/3) and the v_e formula (2)
        EIGHT_THIRDS = TitanFrac(8, 3)

        print("The calculation uses the following fractional components:")
        print(f"v_e^2 ≈ {EIGHT_THIRDS} * {G} * {PI} * {RHO_SHELL} * {R_SQUARED}\n")
        
        # Let's perform the calculation in an order that keeps numbers small for as long as possible.
        # Step 1: Combine constants from the geometric formula
        # (8/3) * pi = (8/3) * (3/1) = 24/3, which simplifies to 8/1. Valid.
        print("Step 1: Combining geometric constants...")
        term1 = EIGHT_THIRDS * PI
        print(f"   {EIGHT_THIRDS} * {PI} = {term1}")
        
        # Step 2: Combine density and radius terms
        # rho_shell * R^2 = (3/1 * 10^2) * (4/1 * 10^12) = 12/1 * 10^14. Valid.
        print("Step 2: Combining density and radius terms...")
        term2 = RHO_SHELL * R_SQUARED
        print(f"   {RHO_SHELL} * {R_SQUARED} = {term2}")

        # Step 3: Multiply the results together with G. This is where the failure occurs.
        # (term1 * G) = (8/1) * (13/2 * 10^-11)
        # This operation results in (104/2), simplifying to (52/1).
        # The numerator 52 is greater than 15, violating the constraint.
        print("Step 3: Multiplying with the gravitational constant G...")
        term1_G = term1 * G
        # The script will fail on the line above and jump to the except block.
        
    except TitanError as e:
        print("-" * 50)
        print("SIMULATION FAILED.")
        print(f"Reason: {e}")
        print("\nCONCLUSION:")
        print("The calculation is not possible. The magnitude of the physical constants (G) and planetary dimensions (R, rho), even with favorable approximations, leads to intermediate products that exceed the 4-bit representation limit (0-15).")
        print("The architecture's constraints prevent the handling of numbers like '52', which arise unavoidably during the calculation.")

solve()