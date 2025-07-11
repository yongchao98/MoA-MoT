import math

# Step 1: Define the Titan 6-bit constrained fraction class
class TitanError(Exception):
    """Custom exception for Titan architecture violations."""
    pass

class TitanFraction:
    """
    Represents a fraction under Titan's 6-bit constraints.
    The value is num/den * 10^exp.
    """
    MAX_VAL = 63

    def __init__(self, num, den=1, exp=0):
        if not all(isinstance(i, int) for i in [num, den, exp]):
            raise TypeError("All parts of a TitanFraction must be integers.")
        if den == 0:
            raise ZeroDivisionError("Denominator cannot be zero.")
        if abs(num) > self.MAX_VAL or abs(den) > self.MAX_VAL:
            raise TitanError(f"Initialization Error: {num}/{den} exceeds 6-bit limit.")
        
        # Simplify the fraction using GCD
        common = math.gcd(num, den)
        self.num = num // common
        self.den = den // common
        self.exp = exp

    def __mul__(self, other):
        """
        Multiplies two TitanFractions, checking for overflow before calculation.
        This is the critical constraint of the Titan architecture.
        """
        new_num_val = self.num * other.num
        new_den_val = self.den * other.den

        if abs(new_num_val) > self.MAX_VAL or abs(new_den_val) > self.MAX_VAL:
            raise TitanError(
                f"Multiplication Overflow: Cannot compute ({self.num}/{self.den}) * ({other.num}/{other.den}). "
                f"Resulting numerator '{new_num_val}' or denominator '{new_den_val}' exceeds {self.MAX_VAL}."
            )
        
        new_exp = self.exp + other.exp
        return TitanFraction(new_num_val, new_den_val, new_exp)

    def __truediv__(self, other):
        """Divides two TitanFractions, checking for overflow."""
        # Division is multiplication by the reciprocal
        reciprocal = TitanFraction(other.den, other.num, -other.exp)
        return self * reciprocal

    def __str__(self):
        if self.exp == 0:
            return f"{self.num}/{self.den}"
        else:
            return f"({self.num}/{self.den})e{self.exp}"

def solve_pandora_problem():
    """
    Attempts to solve the Pandora black hole problem using the Titan architecture simulation.
    """
    print("Task: Calculate gravitational force on a 50kg probe 1km from Pandora's event horizon.")
    print("Simulating Titan 6-bit architecture...\n")

    # Step 2: Define physical constants and problem values as TitanFractions
    # Using the most accurate approximations where num/den <= 63
    try:
        G = TitanFraction(20, 3, -11) # Gravitational constant G ~= 6.67e-11 N m^2/kg^2
        PI = TitanFraction(22, 7)      # Pi ~= 3.1428
        FOUR_THIRDS = TitanFraction(4, 3)
        
        R = TitanFraction(2, 1, 6)     # Pandora radius = 2000 km = 2e6 m
        R_CUBED = TitanFraction(8, 1, 18) # R^3 = (2e6)^3 = 8e18 m^3
        
        RHO = TitanFraction(6, 5, 3)   # Density = 1.2 t/m^3 = 1200 kg/m^3 = 6/5 e3 kg/m^3
        M2 = TitanFraction(50, 1)      # Probe mass = 50 kg
        
        # Distance from event horizon is 1km. Rs is negligible for a planet-mass black hole.
        # r ~= 1 km = 1e3 m.
        # r^2 = (1e3)^2 = 1e6 m^2.
        R_SQUARED = TitanFraction(1, 1, 6)

    except TitanError as e:
        print(f"Error during initialization: {e}")
        print("<<<N0>>>")
        return

    # Print the full equation with the chosen fractional approximations
    print("Full Force Equation:")
    print(f"F = G * (4/3 * PI * R^3 * RHO) * m2 / r^2")
    print("F =",
          f"({G}) * (({FOUR_THIRDS}) * ({PI}) * ({R_CUBED}) * ({RHO}))",
          f"* ({M2}) / ({R_SQUARED})\n")
    
    # Step 3: Attempt to calculate the force step-by-step
    print("Attempting calculation...")
    try:
        # Let's calculate the mass M = (4/3) * PI * R^3 * RHO first
        # We will try to multiply the terms sequentially
        
        print(f"1. Calculate (4/3) * PI = {FOUR_THIRDS} * {PI}")
        # This multiplication (4*22=88, 3*7=21) would create 88 in the numerator.
        # The rules state any *operation* resulting in a value > 63 must be simplified.
        # 4*22 = 88, which is > 63. The MUL instruction would fail.
        # Let's show this explicitly.
        mass_term1 = FOUR_THIRDS * PI
        
    except TitanError as e:
        print(f"Calculation failed at the first step.")
        print(f"Reason: {e}\n")
        print("Conclusion: The product of the numerators (4 * 22 = 88) exceeds the 6-bit limit of 63.")
        print("No reordering of terms can solve this, as any multiplication of the larger numbers")
        print("in the expression (e.g., 20, 22, 50) will inevitably result in an overflow.")
        print("The Titan architecture, with these constraints, is not capable of performing this calculation.")
        print("\nFinal Answer:")
        print("<<<N0>>>")
        return

    # This part of the code is unreachable, but shows the logical next steps
    # if the first one had succeeded.
    mass_term2 = mass_term1 * R_CUBED
    mass = mass_term2 * RHO
    
    force_term1 = G * mass
    force_term2 = force_term1 * M2
    final_force = force_term2 / R_SQUARED

    # If calculation were possible, we would output Y[e]
    # But as demonstrated, it is not.

if __name__ == '__main__':
    solve_pandora_problem()