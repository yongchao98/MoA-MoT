import sys

def display_fermionic_partition_function():
    """
    This script prints the formula for the fermionic partition function Z
    in the imaginary time representation using Feynman’s path integral formalism.
    """
    # Print the main formula
    print("The formula for the fermionic partition function (Z) is given by the path integral over Grassmann fields:")
    print("\nZ = \u222B D\u03C8\u0304 D\u03C8 exp(-S[\u03C8\u0304, \u03C8])\n")

    # Print the formula for the Action S
    print("Where the Euclidean action (S) is:")
    print("\nS[\u03C8\u0304, \u03C8] = \u222B\u2080^\u03B2 d\u03C4 \u222B d\u207Dx [ \u03C8\u0304(\u2202_\u03C4 + H\u2080 - \u03BC)\u03C8 + L_int(\u03C8\u0304, \u03C8) ]\n")

    # Explain the components as requested
    print("--- Components of the Final Equations ---")
    
    print("\nIn the equation Z = \u222B D\u03C8\u0304 D\u03C8 exp(-S[\u03C8\u0304, \u03C8]):")
    print("Z: The partition function of the fermionic system.")
    print("=: Equals sign.")
    print("\u222B: The functional integral symbol, representing integration over all possible field configurations.")
    print("D\u03C8\u0304 D\u03C8: The integration measure for the path integral over the Grassmann fields \u03C8\u0304 and \u03C8.")
    print("exp(...): The exponential function.")
    print("-: The negative sign, indicating e to the power of negative S.")
    print("S[\u03C8\u0304, \u03C8]: The Euclidean action, which is a functional of the fields.")

    print("\nIn the action S = \u222B d\u03C4 \u222B d\u207Dx [ ... ]:")
    print("S: The Euclidean action.")
    print("\u222B\u2080^\u03B2 d\u03C4: Integral over imaginary time \u03C4 from 0 to \u03B2.")
    print("\u222B d\u207Dx: Integral over d-dimensional spatial coordinates x.")
    print("\u03C8\u0304 and \u03C8: The anticommuting Grassmann fields that represent the fermions.")
    print("\t- These fields must satisfy anti-periodic boundary conditions in imaginary time: \u03C8(x, \u03C4) = -\u03C8(x, \u03C4 + \u03B2).")
    print("\u2202_\u03C4: The partial derivative with respect to imaginary time \u03C4.")
    print("H\u2080: The single-particle Hamiltonian (e.g., the kinetic energy operator, -\u210F\u00B2\u2207\u00B2/2m).")
    print("\u03BC: The chemical potential.")
    print("L_int: The Lagrangian density describing interactions between fermions (e.g., four-fermion interaction term).")

if __name__ == '__main__':
    # On some systems, especially Windows, the default console might not support Unicode characters.
    # This attempts to set the encoding to UTF-8 for better character display.
    if sys.stdout.encoding != 'utf-8':
        try:
            sys.stdout.reconfigure(encoding='utf-8')
        except TypeError:
            # This can happen in some environments like IDLE.
            # A fallback will just print without Unicode.
            print("Warning: Could not reconfigure stdout to utf-8. Mathematical symbols may not display correctly.", file=sys.stderr)
            
    display_fermionic_partition_function()
    
    final_answer = "Z = \u222B D\u03C8\u0304 D\u03C8 exp(-\u222B\u2080^\u03B2 d\u03C4 \u222B d\u207Dx [ \u03C8\u0304(\u2202_\u03C4 + H\u2080 - \u03BC)\u03C8 + L_int ])"
    sys.stdout.flush()
    print(f"\n<<<{final_answer}>>>")