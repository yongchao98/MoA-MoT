import sys

def display_fermionic_partition_function_formula():
    """
    This function prints the formula for the fermionic partition function Z
    in the imaginary time path integral formalism.
    """
    print("The formula for the fermionic partition function (Z) in the imaginary time representation is:")
    print("========================================================================================\n")
    
    # Main formula
    print("Z = Integral[over A.P.B.C.] D(psi_bar) D(psi) * exp(-S)\n")
    
    print("Where:")
    print("  - Z: The partition function.")
    print("  - D(psi_bar) D(psi): The functional integral over all configurations of the fermionic (Grassmann) fields.")
    print("  - S: The Euclidean action of the system.")
    print("  - A.P.B.C.: Stands for Anti-Periodic Boundary Conditions, which are crucial for fermions.\n")

    # Euclidean Action (S)
    print("The Euclidean action S is the integral of the Euclidean Lagrangian (L_E) over imaginary time and space:")
    # Using explicit numbers 0 and beta for the integral limits as requested
    print("S = Integral[from tau=0 to beta] d(tau) * Integral[over all space] d^3(x) * L_E\n")
    
    print("  - tau: Imaginary time.")
    print("  - beta: The inverse temperature, defined as 1 / (k_B * T).")
    print("  - 0: The starting point of the imaginary time integral.")
    print("  - beta: The endpoint of the imaginary time integral.")
    print("  - L_E: The Euclidean Lagrangian density.\n")
    
    # Euclidean Lagrangian (L_E)
    print("The Euclidean Lagrangian density L_E for a general fermionic system is:")
    print("L_E = psi_bar(x, tau) * (d/d(tau) + H) * psi(x, tau)\n")
    
    print("  - psi(x, tau), psi_bar(x, tau): The anti-commuting Grassmann fields.")
    print("  - d/d(tau): The derivative with respect to imaginary time.")
    print("  - H: The Hamiltonian operator for the system (e.g., kinetic + potential energy terms).\n")

    # Anti-Periodic Boundary Conditions
    print("The Anti-Periodic Boundary Conditions (A.P.B.C.) for the fermion fields are:")
    print("psi(x, tau=0) = -psi(x, tau=beta)")

if __name__ == '__main__':
    display_fermionic_partition_function_formula()
