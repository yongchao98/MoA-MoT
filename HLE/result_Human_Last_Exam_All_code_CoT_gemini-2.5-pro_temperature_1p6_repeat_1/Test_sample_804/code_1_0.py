class TitanComputerError(Exception):
    """Custom exception for Titan architecture constraint violations."""
    pass

class Fraction:
    """
    A class to represent and operate on fractions under Titan's 5-bit constraints.
    """
    MAX_INT = 31

    def __init__(self, numerator, denominator):
        if not (0 <= numerator <= self.MAX_INT and 0 <= denominator <= self.MAX_INT):
            raise TitanComputerError(f"Initialization failed: Numerator ({numerator}) or denominator ({denominator}) exceeds 5-bit limit of {self.MAX_INT}.")
        if denominator == 0:
            raise TitanComputerError("Initialization failed: Denominator cannot be zero.")
            
        self.num = numerator
        self.den = denominator

    def __mul__(self, other):
        """
        Multiplies two fractions, checking for 5-bit overflow before returning result.
        """
        new_num = self.num * other.num
        new_den = self.den * other.den

        print(f"Attempting multiplication: ({self.num}/{self.den}) * ({other.num}/{other.den})")
        
        # This is the crucial check for the Titan architecture constraint
        if new_num > self.MAX_INT or new_den > self.MAX_INT:
            print(f"Resulting fraction would be: {new_num}/{new_den}")
            raise TitanComputerError(f"Calculation failed: Resulting numerator ({new_num}) or denominator ({new_den}) exceeds 5-bit limit of {self.MAX_INT}.")
        
        return Fraction(new_num, new_den)

    def __str__(self):
        return f"{self.num}/{self.den}"

def calculate_pandora_force():
    """
    Attempts to calculate Pandora's gravity, demonstrating the Titan computer's limitations.
    """
    print("Evaluating feasibility of calculating Pandora's gravitational force on Titan computer.")
    print("-" * 20)
    print("The calculation requires finding the planet's mass, which depends on the volumes of its core and shell.")
    print("The volume of a sphere is proportional to its radius cubed (R^3).")
    print("To simplify, we can use ratios of the physical parameters.")

    # Define the ratio of the core radius to the outer radius
    # R_core = 50 km, R_outer = 1000 km. Ratio = 50/1000 = 1/20.
    # This ratio is representable in the 5-bit system.
    try:
        r_ratio_str = "1/20"
        print(f"\nThe ratio of the core radius to the outer radius is {r_ratio_str}.")
        r_core_over_router = Fraction(1, 20)
        print(f"This can be represented as the 5-bit fraction: {r_core_over_router}")

        # The mass calculation involves cubing this ratio. We must do this with fractional arithmetic.
        # Let's attempt the first step: squaring the ratio.
        print("\nTo find the volume ratio, we must cube this fraction.")
        print("First step is squaring it: (1/20) * (1/20)")
        
        # This multiplication will fail because 20 * 20 = 400, which is > 31.
        r_ratio_squared = r_core_over_router * r_core_over_router

    except TitanComputerError as e:
        print("\n" + "#" * 12 + " ERROR " + "#" * 12)
        print(e)
        print("#" * 31)
        print("\nThe intermediate calculation to find the volume ratio resulted in a denominator (400) that cannot be stored in a 5-bit register.")
        print("Approximation is not viable here, as the next smallest denominator would introduce massive error.")
        print("Conclusion: The calculation is not feasible with the Titan architecture due to its limited numerical range.")
        print("\nFinal Answer: N0")


# Run the simulation
if __name__ == "__main__":
    calculate_pandora_force()
<<<N0>>>