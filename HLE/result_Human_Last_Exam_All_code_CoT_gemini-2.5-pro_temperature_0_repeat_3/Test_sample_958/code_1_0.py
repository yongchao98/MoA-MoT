def solve_evanescent_energy():
    """
    This function presents the final formulas for the time-averaged stored energy
    per unit area for the electric and magnetic fields of an evanescent wave,
    based on the provided answer choices.
    """
    print("The selected answer is D. The formulas are as follows:")
    print("=" * 70)

    # --- Electric Field Energy ---
    print("Final equation for Energy in E field:")
    print("Numerator: n^2 * (2*n^2*sin(theta)^2 - 1)")
    print("Denominator: 2 * (omega/c) * (n^2 - 1) * ((n^2 + 1)*sin(theta)^2 - 1) * sqrt(n^2*sin(theta)^2 - 1)")
    print("Factors: epsilon_0 * |E_x0_i|^2")
    
    print("\nNumbers present in the E-field equation:")
    print("From numerator (2*n^2...): 2")
    print("From numerator (... - 1): 1")
    print("From denominator start (2 * ...): 2")
    print("From denominator (n^2 - 1): 1")
    print("From denominator ((n^2 + 1)...): 1")
    print("From denominator (...*sin(theta)^2 - 1): 1")
    print("From denominator sqrt(... - 1): 1")
    
    print("\n" + "-"*40 + "\n")

    # --- Magnetic Field Energy ---
    print("Final equation for Energy in H field:")
    print("Numerator: n^2 * (n^2*sin(theta)^2 - 1)")
    print("Denominator: 2 * (omega/c) * (n^2 - 1) * ((n^2 + 1)*sin(theta)^2 - 1) * sqrt(n^2*sin(theta)^2 - 1)")
    print("Factors: epsilon_0 * |E_x0_i|^2")

    print("\nNumbers present in the H-field equation:")
    print("From numerator (n^2*sin...): (implicitly 1)")
    print("From numerator (... - 1): 1")
    print("From denominator start (2 * ...): 2")
    print("From denominator (n^2 - 1): 1")
    print("From denominator ((n^2 + 1)...): 1")
    print("From denominator (...*sin(theta)^2 - 1): 1")
    print("From denominator sqrt(... - 1): 1")
    
    print("=" * 70)

if __name__ == "__main__":
    solve_evanescent_energy()