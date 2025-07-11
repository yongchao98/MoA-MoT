import sympy

def solve_electron_interaction():
    """
    This function derives and displays the effective electron-electron interaction
    by integrating out phonon fields, as per the user's request.
    """

    # Define the symbols used in the derivation for symbolic representation.
    g, m, q_j, w_q, Omega = sympy.symbols('g m q_j w_q Omega', real=True, positive=True)

    # --- Derivation Steps ---
    # The effective interaction potential V_eff(q, j, Omega) is derived from the
    # path integral over the phonon fields. The result from a second-order
    # cumulant expansion of the action is:
    #
    # V_eff = (g_qj * g_{-q,j}) * D(q, Omega)
    #
    # where g_qj is the vertex coupling factor and D(q, Omega) is the phonon propagator.

    # 1. The coupling factor product g_qj * g_{-q,j}
    # From the interaction term g * (i*q_j) * ..., we have:
    # g_qj is proportional to i*q_j
    # g_{-q,j} is proportional to i*(-q_j)
    # Their product: (i*q_j) * i*(-q_j) = -i^2 * q_j^2 = q_j^2.
    # Including all constants from the Hamiltonian:
    coupling_product_str = f"g**2 * q_j**2 / (2 * m * w_q)"
    coupling_product = g**2 * q_j**2 / (2 * m * w_q)

    # 2. The phonon propagator D(q, Omega)
    # The propagator for the phonon field phi_q = a_q + a_{-q}^dagger is:
    # D(q, Omega) = 2*w_q / (w_q**2 + Omega**2)
    propagator_str = f"2 * w_q / (w_q**2 + Omega**2)"
    propagator = 2*w_q / (w_q**2 + Omega**2)

    # 3. The effective interaction potential V_eff
    # V_eff is the product of the above terms.
    V_eff = sympy.simplify(coupling_product * propagator)

    # --- Final Output ---
    print("The effective electron-electron interaction potential V_eff for a mode (q, j) is:")
    print("-" * 70)

    # To satisfy the user's request to see each number, we build the equation string explicitly.
    final_eq_str = f"V_eff(q, j, Omega) = (g**2 * q_j**2 / m) * (1 / (w_q**2 + Omega**2))"
    print(final_eq_str)
    print("-" * 70)

    print("\nHere is a breakdown of the equation's components:")
    print(f"\n1. Contribution from the squared electron-phonon coupling vertices:")
    print(f"   (g * i*q_j / (2*m*w_q)**(1/2)) * (g * i*(-q_j) / (2*m*w_q)**(1/2)) * (electron densities)")
    print(f"   This simplifies to a factor of: g**2 * q_j**2 / (2 * m * w_q)")

    print(f"\n2. Contribution from the exchanged phonon's dynamics (propagator):")
    print(f"   The propagator D(q, Omega) for the phonon field is: {propagator_str}")

    print(f"\n3. Combining these parts gives the final potential V_eff:")
    print(f"   V_eff = ( {coupling_product_str} ) * ( {propagator_str} )")
    print(f"   V_eff = {V_eff}")

    print("\nNote: The interaction is repulsive because the imaginary coupling constant squared is negative, i.e., (i)^2 = -1, which cancels the minus sign that typically leads to attraction in phonon-mediated interactions.")

if __name__ == '__main__':
    solve_electron_interaction()