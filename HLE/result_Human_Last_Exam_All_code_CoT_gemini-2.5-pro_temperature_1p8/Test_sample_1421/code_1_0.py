import math

def display_fermionic_partition_function():
    """
    This function prints the formula for the fermionic partition function Z
    in the imaginary time path integral formalism.
    """
    print("The formula for the fermionic partition function Z in the imaginary time representation is a path integral over Grassmann fields.")
    print("=" * 90)

    # Part 1: The Partition Function Integral
    print("\n1. The Partition Function (Z)\n")
    print("   Z is expressed as a functional integral over all possible configurations of the fermionic")
    print("   Grassmann fields psi_bar and psi, with each configuration weighted by e^(-S_E).\n")
    print("      Z = INTEGRAL [D(psi_bar) D(psi)] * exp(-S_E[psi_bar, psi])")
    print("\n   - [D(psi_bar) D(psi)] represents the functional integration measure.")
    print("   - S_E is the Euclidean action of the system.")
    
    print("\n" + "-" * 90 + "\n")

    # Part 2: The Euclidean Action
    print("2. The Euclidean Action (S_E)\n")
    print("   The action describes the system's dynamics in imaginary time (tau). A general form is:")
    print("\n      S_E = INTEGRAL from 0 to beta [d(tau)] * INTEGRAL [d^d(x)] * L_E")
    print("\n   Where the Euclidean Lagrangian density, L_E, for a system in the grand canonical ensemble is often:")
    print("\n      L_E = psi_bar(x, tau) * [ (d/d(tau)) + H_0 - mu ] * psi(x, tau)")
    print("\n   The components are:")
    print(f"     - Integration over imaginary time 'tau' is from 0 to beta.")
    print(f"     - beta = 1 / (k_B * T), where T is temperature and k_B is the Boltzmann constant.")
    print(f"     - H_0 is the single-particle Hamiltonian operator (e.g., kinetic energy).")
    print(f"     - mu is the chemical potential.")

    print("\n" + "-" * 90 + "\n")

    # Part 3: Anti-Periodic Boundary Conditions
    print("3. Anti-Periodic Boundary Conditions\n")
    print("   Due to the trace in the definition of Z and fermionic statistics, the fields are")
    print("   anti-periodic in imaginary time.\n")
    print("      psi(x, beta) = -psi(x, 0)")
    print("      psi_bar(x, beta) = -psi_bar(x, 0)")
    print("\n   This is a fundamental difference compared to bosonic systems, which have periodic boundary conditions.")
    print("\n" + "=" * 90)

if __name__ == "__main__":
    display_fermionic_partition_function()