import numpy as np

def calculate_effective_interaction(g, m, w_q, q_vector):
    """
    Calculates the effective electron-electron interaction potential mediated by phonons.

    Args:
        g (float): The electron-phonon coupling constant.
        m (float): The mass parameter from the coupling term.
        w_q (float): The phonon frequency for the given momentum q.
        q_vector (list or np.ndarray): The phonon momentum vector q.
    """
    q_vector = np.array(q_vector)
    
    # Calculate the squared magnitude of the momentum vector q
    q_squared = np.dot(q_vector, q_vector)

    # Calculate the effective interaction potential V_eff(q)
    # V_eff(q) = - (g^2 * |q|^2) / (m * w_q^2)
    v_eff = - (g**2 * q_squared) / (m * w_q**2)

    # Print the step-by-step calculation and the final result
    print("Derivation of the effective electron-electron interaction potential V_eff(q):")
    print("-----------------------------------------------------------------------")
    print(f"The general formula is: V_eff(q) * rho(q) * rho(-q), where V_eff(q) = - (g^2 * |q|^2) / (m * w_q^2)")
    print("\nGiven parameters:")
    print(f"  g (coupling constant) = {g}")
    print(f"  m (mass parameter) = {m}")
    print(f"  w_q (phonon frequency) = {w_q}")
    print(f"  q (momentum vector) = {q_vector}")
    
    print("\nCalculation steps:")
    print(f"  1. Calculate |q|^2:")
    print(f"     |q|^2 = {' + '.join([f'{qi}^2' for qi in q_vector])} = {' + '.join([str(qi**2) for qi in q_vector])} = {q_squared}")
    
    print(f"\n  2. Substitute the values into the V_eff(q) formula:")
    # Numerator and denominator calculation
    numerator = g**2 * q_squared
    denominator = m * w_q**2
    print(f"     V_eff(q) = - ({g}^2 * {q_squared}) / ({m} * {w_q}^2)")
    print(f"     V_eff(q) = - ({numerator}) / ({denominator})")
    print(f"     V_eff(q) = {v_eff}")

    print("\nFinal effective interaction Hamiltonian term for the given q:")
    print(f"  H_eff(q) = ({v_eff}) * rho(q) * rho(-q)")


# --- Example Usage ---
# Define the parameters for the calculation
g_const = 1.5      # Coupling constant
mass_param = 0.5   # Mass from the coupling term
omega_q = 20.0     # Phonon frequency
q_vec = [1.0, 1.0, 2.0] # Momentum vector q

# Calculate and display the result
calculate_effective_interaction(g_const, mass_param, omega_q, q_vec)
