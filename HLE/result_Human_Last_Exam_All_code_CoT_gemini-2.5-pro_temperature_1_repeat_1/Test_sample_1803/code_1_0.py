import sympy

def solve_energy_shift():
    """
    Calculates the ground state energy shift for two interacting quantum harmonic oscillators
    using symbolic mathematics.
    """
    # --- Step 1: Define all symbols ---
    # We use sympy for symbolic mathematics.
    # Physical constants (real and positive)
    e, m, w0, hbar = sympy.symbols('e m omega_0 hbar', real=True, positive=True)
    # Coordinates and distance (real)
    x1, x2, R = sympy.symbols('x_1 x_2 R', real=True)

    print("--- Task: Ground State Energy Shift of Two Interacting QHOs ---")
    print("\nThis script calculates the energy shift using the following plan:")
    print("1. Define the Coulomb interaction potential V between the four charges.")
    print("2. Taylor expand V for large R (R >> x1, x2) to find the leading interaction term (V_pert).")
    print("3. Apply perturbation theory to find the ground state energy shift.\n")

    # --- Step 2: Define and Expand the Interaction Potential V ---
    # The problem statement uses the Coulomb interaction e^2/(4*pi*r), and for simplicity,
    # we can treat the constant 1/(4*pi) as 1, leading to an interaction of e^2/r.
    # Positions of charges:
    # Oscillator 1 (+e, -e): located at x1/2 and -x1/2 relative to the origin.
    # Oscillator 2 (+e, -e): located at R + x2/2 and R - x2/2 relative to the origin.
    
    # Distances between the four charges
    r_pp = R + (x2 - x1) / 2  # (+e of osc1) to (+e of osc2)
    r_pm = R - (x2 + x1) / 2  # (+e of osc1) to (-e of osc2)
    r_mp = R + (x2 + x1) / 2  # (-e of osc1) to (+e of osc2)
    r_mm = R - (x2 - x1) / 2  # (-e of osc1) to (-e of osc2)

    # Full Coulomb potential
    V = e**2 * (1/r_pp - 1/r_pm - 1/r_mp + 1/r_mm)

    # To perform the Taylor expansion for R >> x1, x2, we expand in x1 and x2 around 0.
    # We expand up to a sufficiently high order to capture the first non-zero term.
    V_expanded = sympy.series(V, x1, 0, 3)
    V_expanded = sympy.series(V_expanded.removeO(), x2, 0, 3)
    
    # Simplify the resulting expression to get the leading term.
    V_pert = sympy.simplify(V_expanded.removeO())

    print("--- Step 1 & 2: Interaction Potential ---")
    print("The leading term of the interaction potential (the perturbation V_pert) is:")
    # Using pprint for better formatting of the symbolic expression
    sympy.pprint(V_pert, use_unicode=True)
    print("\nThis is the classic dipole-dipole interaction potential.\n")

    # --- Step 3: Perturbation Theory ---
    print("--- Step 3: Perturbation Theory Calculation ---")
    
    # 3a. First-order correction
    print("a) First-Order Correction (ΔE^(1)):")
    print("ΔE^(1) = <0,0| V_pert |0,0>, where |0,0> is the ground state.")
    print("This is proportional to <0|x1|0> * <0|x2|0>. Since the expectation value of position <n|x|n> is 0 for any QHO state, the first-order correction is zero.")
    print("ΔE^(1) = 0\n")

    # 3b. Second-order correction
    print("b) Second-Order Correction (ΔE^(2)):")
    print("The formula is: ΔE^(2) = Σ_{k≠0} |<k| V_pert |0,0>|^2 / (E_00 - E_k)")
    
    # The matrix element <n1,n2| V_pert |0,0> is proportional to <n1|x1|0><n2|x2|0>.
    # The QHO matrix element <n|x|0> is non-zero only for n=1.
    # Therefore, the sum reduces to a single term with the intermediate state k = |1,1>.
    
    # Define the required QHO matrix element <1|x|0>
    x_01 = sympy.sqrt(hbar / (2 * m * w0))
    print("The perturbation only connects the ground state |0,0> to the excited state |1,1>.")
    print("The position matrix element is <1|x|0> = ", end="")
    sympy.pprint(x_01, use_unicode=True)
    
    # Calculate the perturbation matrix element M = <1,1| V_pert |0,0>
    # V_pert = C * x1 * x2, so M = C * <1|x1|0> * <1|x2|0>
    C = V_pert / (x1 * x2)
    M = C * x_01 * x_01
    
    print("\nThe perturbation matrix element M = <1,1| V_pert |0,0> is:")
    sympy.pprint(M, use_unicode=True)

    # The energy denominator is E_00 - E_11
    # E_n = (n + 1/2)*hbar*w0 => E_00 = hbar*w0, E_11 = 3*hbar*w0
    energy_denom = -2 * hbar * w0
    print("\nThe energy denominator (E_00 - E_11) is:")
    sympy.pprint(energy_denom, use_unicode=True)

    # Calculate the final energy shift
    E_shift = (M**2) / energy_denom
    E_shift_simplified = sympy.simplify(E_shift)
    
    print("\n--- Final Result: Ground State Energy Shift ---")
    print("The second-order energy shift ΔE = |M|^2 / (E_00 - E_11).")
    print("Substituting the terms we calculated:")
    
    # Construct the final equation string as requested
    num_str = sympy.srepr(sympy.simplify(M**2))
    den_str = sympy.srepr(energy_denom)
    
    print(f"\nΔE = ({sympy.pretty(sympy.simplify(M**2), use_unicode=True)}) / ({sympy.pretty(energy_denom, use_unicode=True)})")
    
    print("\nWhich simplifies to the final result:")
    sympy.pprint(E_shift_simplified, use_unicode=True)

    return E_shift_simplified

if __name__ == '__main__':
    final_answer = solve_energy_shift()
    # The final answer is returned in the required format
    print(f"\n<<<{final_answer}>>>")