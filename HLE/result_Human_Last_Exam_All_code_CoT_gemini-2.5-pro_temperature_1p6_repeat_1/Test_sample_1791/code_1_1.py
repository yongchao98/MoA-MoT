def solve_electron_interaction():
    """
    This function prints the derived effective electron-electron interaction potential.
    Each part of the equation is printed component by component for clarity.
    """
    print("The effective electron-electron interaction potential V_eff for a given momentum transfer q is derived by integrating out the phonon fields.")
    print("The result is a frequency-dependent potential. Here are the components of the final equation:")
    print("-" * 30)

    # Define the symbols used in the equation
    V_eff = "V_eff(q, i*nu_n)"
    sum_j = "sum over j (phonon polarizations)"
    g = "g (electron-phonon coupling constant)"
    q_j = "q_j (j-th component of momentum q)"
    m = "m (ion mass)"
    nu_n = "nu_n (bosonic Matsubara frequency for energy transfer)"
    w_q = "w_q (phonon frequency at momentum q)"

    # Print each component
    print(f"Interaction Potential: {V_eff}")
    print("This potential is composed of:")
    print(f"  1. A sum over phonon polarizations: {sum_j}")
    print(f"  2. A numerator representing the squared interaction strength for each mode:")
    print(f"     Numerator = {g}^2 * {q_j}^2 / {m}")
    print(f"  3. A denominator from the phonon propagator, dependent on energy and momentum:")
    print(f"     Denominator = {nu_n}^2 + {w_q}^2")
    print(f"  4. An overall negative sign, indicating an attractive interaction.")
    
    print("-" * 30)
    print("The final equation for the interaction potential is:")

    final_equation = f"{V_eff} = - {sum_j} [ ({g}^2 * {q_j}^2 / {m}) / ({nu_n}^2 + {w_q}^2) ]"
    print(final_equation)

if __name__ == '__main__':
    solve_electron_interaction()