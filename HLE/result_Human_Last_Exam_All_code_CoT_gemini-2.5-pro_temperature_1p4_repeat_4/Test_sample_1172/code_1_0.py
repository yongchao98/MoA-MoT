import math

def display_inductance_change_formula():
    """
    This function prints the derived formula for the change
    in mutual inductance per unit length (Δm).
    """

    # Define symbolic representations for the variables
    mu_0 = "μ₀"
    h = "h"
    d = "d"
    R1 = "R₁"
    R2 = "R₂"
    pi = "π"
    delta_m = "Δm"

    # Define the numerical constants in the formula
    numerator_constant = 1  # Implied constant
    denominator_constant = 2
    power_h = 2
    power_d = 2
    power_ratio = 2
    subtracted_constant = 1

    # Print the final formula, highlighting the numerical components as requested.
    print("The change in mutual inductance per unit length is given by the formula:")
    
    final_formula = (
        f"{delta_m} = ({mu_0} * {h}**{power_h} / ({denominator_constant} * {pi} * {d}**{power_d})) * "
        f"(({R2}/{R1})**{power_ratio} - {subtracted_constant})"
    )
    
    print(final_formula)

if __name__ == "__main__":
    display_inductance_change_formula()