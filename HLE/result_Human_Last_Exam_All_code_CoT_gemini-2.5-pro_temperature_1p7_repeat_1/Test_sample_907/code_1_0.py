def print_absorption_equation(case='a'):
    """
    This function prints the equation for the absorption cross-section
    for a chain of molecules based on the specified case.

    Args:
        case (str): 'a' for non-interacting molecules,
                    'b' for nearest-neighbor interaction.
    """
    if case == 'a':
        print("--- Case a) The interaction between molecules can be neglected. ---\n")
        print("In this case, the N molecules are independent. The final states of the transition")
        print("are N degenerate states, where each state |n> corresponds to the n-th molecule")
        print("being excited. The absorption is the incoherent sum of absorptions by each molecule.")
        print("\nThe general equation for the absorption cross-section sigma(omega_L) is:")
        print("sigma(omega_L) = C * Sum_f [ (Ef - Ei) * |<f|d|i>|^2 * exp(-( ((Ef - Ei)/hbar - omega_L)^2 * tau^2 )) ]\n")
        print("For non-interacting molecules, this simplifies because the terms in the sum are identical:")
        print("  - Initial state |i> is the ground state |G> with energy Ei = 0.")
        print("  - Final states |f> are N states |n> (n=1...N) with the same energy Ef = epsilon_e.")
        print("  - The transition dipole moment <n|d|G> is mu for each molecule n.")
        print("\nThe final equation is:")
        print("sigma(omega_L) = C_const * N * epsilon_e * |mu|^2 * exp(-( (epsilon_e/hbar - omega_L)^2 * tau^2 ))")
        print("\nWhere the terms are:")
        print("  C_const  : A constant prefactor combining physical constants.")
        print("  N        : The number of molecules in the chain.")
        print("  epsilon_e: The excitation energy of a single molecule.")
        print("  mu       : The transition dipole moment of a single molecule.")
        print("  hbar     : The reduced Planck's constant.")
        print("  omega_L  : The central frequency of the Gaussian laser pulse.")
        print("  tau      : The duration of the Gaussian laser pulse.")
        print("-" * 70)

    elif case == 'b':
        print("--- Case b) The interaction between near-neighbors should be considered. ---\n")
        print("With nearest-neighbor coupling J, the excitations delocalize into a band of Frenkel exciton states.")
        print("The final states are no longer degenerate but form an excitonic band.")
        print("\nThe general equation for the absorption cross-section sigma(omega_L) is a sum over these exciton states |k>:")
        print("sigma(omega_L) = C * Sum_{k=1 to N} [ Ek * |<k|d|G>|^2 * exp(-( (Ek/hbar - omega_L)^2 * tau^2 )) ]\n")
        print("For a chain with open boundary conditions, the terms are:")
        print("  |G>                : The initial ground state with energy E_ground = 0.")
        print("  |k>                : The final exciton state, with k = 1, 2, ..., N.")
        print("  Ek                 : The energy of the exciton state |k>.")
        print("                       Ek = epsilon_e + 2*J*cos(k*pi / (N+1))")
        print("  <k|d|G>            : The transition dipole moment from the ground state to the exciton state |k>.")
        print("                       <k|d|G> = mu * sqrt(2/(N+1)) * Sum_{n=1 to N}[sin(n*k*pi / (N+1))]")
        print("  Selection Rule     : The sum for <k|d|G> is non-zero only for odd values of k. Thus, only")
        print("                       transitions to exciton states with k=1, 3, 5, ... are optically allowed.")
        print("\nThe final equation is a sum over the allowed transitions (k=odd):")
        print("sigma(omega_L) = C_const * Sum_{k=1,3,...}^{N} [ (epsilon_e + 2*J*cos(k*pi/(N+1))) * |M_k|^2 * exp(-( ((epsilon_e + 2*J*cos(k*pi/(N+1)))/hbar - omega_L)^2 * tau^2 )) ]")
        print("with |M_k|^2 = |mu * sqrt(2/(N+1)) * Sum_{n=1 to N}[sin(n*k*pi/(N+1))]|^2")
        print("\nWhere the additional terms are:")
        print("  J        : The nearest-neighbor coupling energy.")
        print("  k        : The exciton quantum number (an integer from 1 to N).")
        print("-" * 70)

    else:
        print("Invalid case specified. Please choose 'a' or 'b'.")


if __name__ == '__main__':
    # Print the equation for case a)
    print_absorption_equation('a')
    # Print the equation for case b)
    print_absorption_equation('b')