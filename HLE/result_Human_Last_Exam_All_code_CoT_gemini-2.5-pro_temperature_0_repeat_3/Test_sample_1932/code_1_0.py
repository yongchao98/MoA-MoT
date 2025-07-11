def display_hourglass_weight_change_formula():
    """
    This function prints the derived formula for the change in weight of a running hourglass.
    The derivation assumes the dominant effect is the impact force of the sand,
    magnified by the Janssen (arching) effect in the top chamber.
    """
    
    # The formula for the change in weight (ΔW) is constructed from:
    # Mass flow rate (ṁ) = (Total Mass) / t = (ρ * π * d² * h / 4) / t
    # Impact velocity (v) = sqrt(2 * g * L)
    # We use the maximum fall distance L = H for the largest possible effect.
    # ΔW ≈ ṁ * v
    
    print("The estimated change in weight, ΔW, considering the largest possible effect, is given by the expression:")
    
    # The prompt requires outputting the numbers/variables in the final equation.
    # We will print the formula in a human-readable text format.
    # The formula corresponds to: (π * d**2 * h * ρ / (4 * t)) * sqrt(2 * g * H)
    
    d = "d"
    h = "h"
    H = "H"
    rho = "ρ"
    t = "t"
    g = "g"
    pi = "π"
    
    print(f"ΔW = ({pi} * {d}**2 * {h} * {rho} / (4 * {t})) * sqrt(2 * {g} * {H})")

display_hourglass_weight_change_formula()