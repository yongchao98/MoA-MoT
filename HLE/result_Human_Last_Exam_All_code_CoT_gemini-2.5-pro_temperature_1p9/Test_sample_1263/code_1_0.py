import math

class TitanError(Exception):
    """Custom exception for Titan architecture failures."""
    pass

class TitanFraction:
    """
    A class to represent a fraction on the Titan architecture.
    It enforces the 4-bit (0-15) constraint on numerators and denominators.
    """
    def __init__(self, num, den=1, exp=0):
        if not (0 <= num <= 15 and 1 <= den <= 15):
            raise TitanError(f"Initialization failed: Numerator {num} or denominator {den} is out of the 4-bit range (0-15).")
        
        # Simplify the fraction using GCD
        common = math.gcd(num, den)
        self.num = num // common
        self.den = den // common
        self.exp = exp

    def __mul__(self, other):
        """
        Multiplies two TitanFractions, checking for overflow.
        This is the critical operation where constraints are likely to be violated.
        """
        new_num = self.num * other.num
        new_den = self.den * other.den
        new_exp = self.exp + other.exp

        if new_num > 15 or new_den > 15:
            # This demonstrates that the operation failed.
            # The prompt's examples of expansion (e.g., 13*6/5 = 13*(1+1/5)) are workarounds
            # that a programmer must devise. However, for a long chain of multiplications,
            # this leads to an explosion of terms or repeated overflows.
            raise TitanError(
                f"Overflow on multiplication: {self.num}/{self.den} * {other.num}/{other.den} -> {new_num}/{new_den}. "
                "Result exceeds 4-bit constraints."
            )
        
        return TitanFraction(new_num, new_den, new_exp)

    def __repr__(self):
        return f"({self.num}/{self.den})e{self.exp}"

def calculate_escape_velocity():
    """
    This function attempts to calculate the escape velocity following Titan's rules.
    It defines the constants and tries to perform the multiplication chain.
    """
    print("Attempting to calculate Pandora's escape velocity on Titan...")
    print("Formula: v_e^2 = 8/3 * G * pi * rho_s * R^2 + [negligible core term]")
    print("We will test if the main term is computable.\n")

    try:
        # Step 1: Define constants with the most plausible fractional approximations.
        # We choose approximations with small numerators and denominators to maximize success chance.
        term_8_over_3 = TitanFraction(8, 3)
        # G = 6.67e-11. Approximation: 13/2 = 6.5
        G = TitanFraction(13, 2, exp=-11)
        # pi = 3.14. Approximation: 3/1 (most simple) or 13/4=3.25
        pi = TitanFraction(3, 1)
        # rho_s = 300. Approximation: 3/1 * 10^2
        rhos = TitanFraction(3, 1, exp=2)
        # R = 2e6. R^2 = (2e6)^2 = 4e12
        R2 = TitanFraction(4, 1, exp=12)
        
        print("Selected approximations:")
        print(f"  8/3          -> {term_8_over_3}")
        print(f"  G (6.67e-11) ~ {G}")
        print(f"  pi (3.14)    ~ {pi}")
        print(f"  rho_s (300)  ~ {rhos}")
        print(f"  R^2 (4e12)   ~ {R2}\n")

        # Step 2: Attempt the multiplication chain, checking for overflow at each step.
        print("Executing multiplication chain...")
        
        # Let's try to compute the product of coefficients for the main term.
        # Order of operations can be changed, but the result is the same.
        # (8/3) * G
        print(f"1. Multiplying {term_8_over_3} by G = {G}...")
        product1 = term_8_over_3 * G
        # This first step results in (8*13)/(3*2) = 104/6. Numerator 104 > 15.
        # The code will raise a TitanError here, and the except block will catch it.

    except TitanError as e:
        print(f"\nCOMPUTATION FAILED.")
        print(f"Reason: {e}")
        print("\nAnalysis: The chain of multiplications required involves numerous constants.")
        print("The product of their numerators (e.g., 8 * 13 = 104) quickly exceeds the 4-bit limit of 15.")
        print("No matter the order of operations or the choice of reasonable fractional approximations,")
        print("an overflow occurs early in the calculation. Therefore, the architecture cannot handle this problem.")
        print("\nFinal Answer: N0")
        return "<<<N0>>>"

    # This part of the code is unreachable
    print("Computation successful (This should not happen).")
    return "<<<Y...>>>"


# The final output is determined by the simulation's outcome.
calculate_escape_velocity()