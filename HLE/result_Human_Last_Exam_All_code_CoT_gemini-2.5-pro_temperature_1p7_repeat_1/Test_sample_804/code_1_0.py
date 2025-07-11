class TitanFraction:
    """
    A class to represent a fraction under Titan's 5-bit constraints.
    The numerator and denominator must be integers between 0 and 31.
    """
    MAX_VAL = 31

    def __init__(self, num, den=1):
        if not (isinstance(num, int) and isinstance(den, int)):
            raise TypeError("Numerator and denominator must be integers.")
        if not (0 <= num <= self.MAX_VAL and 1 <= den <= self.MAX_VAL):
            raise ValueError(
                f"Cannot represent {num}/{den}. Numerator and denominator must be in [0, 31]."
            )
        self.num = num
        self.den = den

    def __repr__(self):
        return f"{self.num}/{self.den}"

def multiply_fractions(f1, f2):
    """
    Multiplies two TitanFractions, checking for overflow before creation.
    This models the constraint that intermediate results must not exceed 31.
    """
    # The result of a multiplication a/b * c/d is (a*c)/(b*d).
    # We must check if the new numerator and denominator are valid.
    new_num = f1.num * f2.num
    new_den = f1.den * f2.den

    if new_num > TitanFraction.MAX_VAL or new_den > TitanFraction.MAX_VAL:
        # Per Titan rules, this operation is illegal as it breaks the 5-bit constraint.
        # We cannot simplify after the fact; the operation itself must be valid.
        # The example in the prompt shows pre-simplification, but 1/20 has no
        # common factors with itself to simplify.
        raise OverflowError(
            f"Multiplying {f1} * {f2} results in {new_num}/{new_den}, which exceeds the 5-bit limit."
        )
    return TitanFraction(new_num, new_den)

def calculate_gravity_on_titan():
    """
    Attempts to calculate Pandora's gravity based on Titan's computational rules.
    """
    print("--- Titan Spacecraft Control System Feasibility Test ---")
    print("Objective: Calculate gravitational force on a 30kg probe 500m above Pandora.")
    print("\n[Step 1] Representing physical quantities as Titan fractions.")

    # The formula for total mass involves the ratio of the radii and densities.
    # M_total = 4/3 * pi * [ (r_core^3 * rho_core) + (r_pandora^3 - r_core^3) * rho_shell ]
    # To manage the large numbers, we can factor out terms:
    # M_total = (4/3*pi*rho_shell*r_pandora^3) * [ (r_core/r_pandora)^3 * (rho_core/rho_shell) + 1 - (r_core/r_pandora)^3 ]
    
    # Let's calculate the dimensionless ratios first.
    # Radius ratio: r_core / r_pandora = 50 km / 1000 km = 1/20
    # Density ratio: rho_core / rho_shell = 1.2 / 0.3 = 4/1
    
    try:
        radius_ratio = TitanFraction(1, 20)
        print(f"Radius ratio (r_core/r_pandora): {radius_ratio}")
        
        density_ratio = TitanFraction(4, 1)
        print(f"Density ratio (rho_core/rho_shell): {density_ratio}")

    except (ValueError, OverflowError) as e:
        print(f"\n[Error] Failed to initialize basic parameters: {e}")
        print("\n[Conclusion] The calculation is NOT feasible.")
        print(">>> N0")
        return

    print("\n[Step 2] Attempting to compute the term (radius_ratio)^2 as part of the mass calculation.")
    print(f"This requires the operation: {radius_ratio} * {radius_ratio}")

    try:
        # This is the critical step that demonstrates the machine's limitation.
        radius_ratio_squared = multiply_fractions(radius_ratio, radius_ratio)
        print(f"Result: {radius_ratio_squared}")
        print("\nThis line should not be reached.")

    except (ValueError, OverflowError) as e:
        print(f"\n[FATAL ERROR] Operation failed. Reason: {e}")
        print("\nThe Titan computer cannot perform this essential calculation.")
        print("To compute the mass, we need to square 1/20. The denominator becomes 20 * 20 = 400.")
        print("Since 400 is greater than 31, the operation violates the 5-bit constraint.")
        print("It's impossible to even represent the physical parameters or perform basic arithmetic with them.")
        print("The final gravitational equation cannot be formed or calculated.")
        
    print("\n[Conclusion] The calculation is NOT feasible on Titan.")
    
# Run the simulation
calculate_gravity_on_titan()