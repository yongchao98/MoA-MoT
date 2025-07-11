import fractions

def calculate_relativistic_shift():
    """
    Calculates the first-order relativistic energy shift for a hydrogen atom.
    
    The user is prompted to enter the principal (n) and angular momentum (l)
    quantum numbers. The script then calculates the energy shift based on the
    standard first-order perturbation theory formula and prints the steps.
    """
    # Given quantum numbers
    n = 3
    l = 2

    print(f"Calculating the first-order relativistic energy shift for n={n}, l={l}.")
    print("Formula: ΔE = - (E_n^2 / (2*m*c^2)) * [4*n / (l + 1/2) - 3]\n")

    # Part 1: The bracketed term
    print("Step 1: Evaluate the term in brackets [4*n / (l + 1/2) - 3]")
    part1_num = 4 * n
    part1_den = l + 0.5
    print(f"  = [4*{n} / ({l} + 0.5) - 3]")
    print(f"  = [{part1_num} / {part1_den} - 3]")
    part1_div = part1_num / part1_den
    print(f"  = [{part1_div:.3f} - 3]")
    part1_result = part1_div - 3
    part1_frac = fractions.Fraction(part1_result).limit_denominator()
    print(f"  = {part1_result:.3f} = {part1_frac.numerator}/{part1_frac.denominator}\n")

    # Part 2: The pre-factor
    print("Step 2: Evaluate the pre-factor E_n^2 / (2*m*c^2) in terms of fundamental constants m, c, α")
    print("  The unperturbed energy E_n is given by: E_n = -m*c^2*α^2 / (2*n^2)")
    pre_factor_den_1 = 2 * n**2
    print(f"  For n={n}, E_{n} = -m*c^2*α^2 / (2*{n}^2) = -m*c^2*α^2 / {pre_factor_den_1}")
    print(f"  Then, E_{n}^2 = (m*c^2*α^2)^2 / ({pre_factor_den_1})^2 = m^2*c^4*α^4 / {pre_factor_den_1**2}")
    pre_factor_den_2 = 2 * (pre_factor_den_1**2)
    print(f"  So, E_{n}^2 / (2*m*c^2) = m*c^2*α^4 / (2 * {pre_factor_den_1**2}) = m*c^2*α^4 / {pre_factor_den_2}\n")

    # Part 3: Combine and simplify
    print("Step 3: Combine the terms and simplify the final coefficient")
    print(f"  ΔE = - (m*c^2*α^4 / {pre_factor_den_2}) * ({part1_frac.numerator}/{part1_frac.denominator})")
    final_num = part1_frac.numerator
    final_den = pre_factor_den_2 * part1_frac.denominator
    print(f"  ΔE = - ({final_num} / {final_den}) * m*c^2*α^4")
    final_frac = fractions.Fraction(final_num, final_den).limit_denominator()
    print(f"  Simplified: ΔE = -({final_frac.numerator}/{final_frac.denominator}) * m*c^2*α^4\n")

    # Final answer print
    print("The final result for the energy shift is:")
    print(f"ΔE = -({final_frac.numerator}/{final_frac.denominator}) * m*c^2*α^4")

if __name__ == "__main__":
    calculate_relativistic_shift()