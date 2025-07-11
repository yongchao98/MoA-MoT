import sys

# Set default encoding to UTF-8 to ensure mathematical symbols print correctly
if sys.stdout.encoding != 'UTF-8':
    sys.stdout.reconfigure(encoding='utf-8')

def display_fermionic_partition_function():
    """
    Prints the formula for the fermionic partition function Z in the
    imaginary time path integral formalism and explains its components.
    """
    
    # Using Unicode characters for mathematical symbols
    integral = "\u222B"
    partial = "\u2202"
    psi_bar = "\u03C8\u0304"  # Represents psi with an overbar
    psi = "\u03C8"
    tau = "\u03C4"
    beta = "\u03B2"
    d_dim = "d"

    print("The formula for the fermionic partition function Z using the imaginary time path integral is:")
    print("-" * 80)
    
    # The full formula is constructed as a string.
    # The term inside the exponential is the Euclidean Action, S_E.
    equation = (
        f"  Z = {integral} D{psi_bar} D{psi} exp[ "
        f"-{integral}_0^{beta} d{tau} "
        f"{integral} d^{d_dim}x  "
        f"{psi_bar}({tau}, x)({partial}_{tau} + H_0){psi}({tau}, x) ]"
    )

    print(equation)
    print("-" * 80)
    print("Explanation of the terms:")
    
    print(f"\n  Z: The Grand Canonical Partition Function.")
    
    print(f"\n  {integral} D{psi_bar} D{psi}: This is the path integral, representing a sum over all possible configurations "
          f"of the fermionic fields.")
          
    print(f"\n  {psi}({tau}, x), {psi_bar}({tau}, x): These are anticommuting Grassmann fields. They are functions of imaginary "
          f"time '{tau}' and spatial position 'x'.")

    print(f"\n  exp[...]: The exponential of the negative Euclidean action for the fermions.")
          
    print(f"\n  {beta}: Represents the inverse temperature, where {beta} = 1 / (k_B * T). The integration over "
          f"imaginary time {tau} runs from 0 to {beta}.")

    print(f"\n  {integral} d^{d_dim}x: Represents the integral over all {d_dim} spatial dimensions.")

    print(f"\n  {partial}_{tau}: The partial derivative with respect to imaginary time {tau}.")

    print(f"\n  H_0: The single-particle Hamiltonian, which typically includes the kinetic energy and "
          f"chemical potential (e.g., H_0 = -∇²/2m - μ).")
          
    print("\n  An essential component of this formula is the anti-periodic boundary condition for fermions "
          f"in the imaginary time direction:")
    print(f"      {psi}(0, x) = -{psi}({beta}, x)")
    
if __name__ == '__main__':
    display_fermionic_partition_function()