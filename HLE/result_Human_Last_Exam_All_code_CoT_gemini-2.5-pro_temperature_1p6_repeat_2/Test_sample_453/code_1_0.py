from fractions import Fraction

def calculate_energy_shift():
    """
    Calculates the second-order relativistic energy shift for a hydrogen atom.
    
    The calculation uses a known closed-form expression for the energy shift,
    which depends on the principal quantum number 'n' and the angular momentum
    quantum number 'l'.
    """
    # Given quantum numbers
    n = 3
    l = 2

    # --- Step 1: Calculate S1 and S2 helper sums ---
    # S1(n, l) = 0.5 * (n^2 + l*(l+1) - 1 - 1/(2*l+1))
    s1 = 0.5 * (n**2 + l * (l + 1) - 1 - 1 / (2 * l + 1))
    
    # S2(n, l) = 0.75*n^2 + 1.5*l*(l+1) - 1.5 - (n^2 + l*(l+1))/(2*(2*l+1))
    s2 = 0.75 * n**2 + 1.5 * l * (l + 1) - 1.5 - (n**2 + l * (l + 1)) / (2 * (2 * l + 1))
    
    # --- Step 2: Calculate the term in the main bracket ---
    # B = 0.5 * S2 - S1 + n/(2*(l+0.5)) - 3/8
    bracket_term = 0.5 * s2 - s1 + n / (2 * (l + 0.5)) - (3.0 / 8.0)
    
    # --- Step 3: Calculate the pre-factor and the final coefficient ---
    # Pre-factor = -1 / (4 * n^4)
    pre_factor = -1.0 / (4 * n**4)
    
    # Final coefficient for mc^2 * alpha^8
    final_coefficient = pre_factor * bracket_term
    
    # --- Step 4: Print the detailed calculation and final result ---
    print("Calculation of the Second-Order Relativistic Energy Shift\n")
    print(f"For a hydrogen atom state with n = {n}, l = {l}:\n")
    
    print(f"1. Calculate helper sums S1 and S2:")
    print(f"   S1({n}, {l}) = 0.5 * ({n}^2 + {l}*({l}+1) - 1 - 1/(2*{l}+1)) = {s1}")
    print(f"   S2({n}, {l}) = 0.75*{n}^2 + 1.5*{l}*({l}+1) - 1.5 - ({n}^2 + {l}*({l}+1))/(2*(2*{l}+1)) = {s2}\n")

    print(f"2. Calculate the bracketed term in the energy shift formula:")
    print(f"   Term = [ 0.5 * S2 - S1 + n/(2*(l+0.5)) - 3/8 ]")
    print(f"        = [ 0.5 * {s2:.4f} - {s1:.4f} + {n}/(2*({l}+0.5)) - 3/8 ]")
    print(f"        = [ {0.5 * s2:.4f} - {s1:.4f} + {n / (2 * (l + 0.5)):.4f} - {3.0/8.0:.4f} ]")
    print(f"        = {bracket_term:.4f}\n")
    
    print(f"3. The full energy shift formula is:")
    print(f"   ΔE = - (m*c^2 * α^8) / (4 * n^4) * (Term)")
    print(f"   ΔE = - (m*c^2 * α^8) / (4 * {n}^4) * ({bracket_term:.4f})")
    print(f"   ΔE = - (m*c^2 * α^8) / ({4 * n**4}) * ({bracket_term:.4f})")
    print(f"   ΔE = ({pre_factor:.6f}) * ({bracket_term:.4f}) * (m*c^2 * α^8)")
    print(f"   ΔE = {final_coefficient:.9f} * (m*c^2 * α^8)\n")

    # Express the final result as a fraction
    final_fraction = Fraction(final_coefficient).limit_denominator()
    print("4. Final answer expressed as a fraction:")
    print(f"   ΔE = {final_fraction.numerator}/{final_fraction.denominator} * m*c^2 * α^8")

calculate_energy_shift()