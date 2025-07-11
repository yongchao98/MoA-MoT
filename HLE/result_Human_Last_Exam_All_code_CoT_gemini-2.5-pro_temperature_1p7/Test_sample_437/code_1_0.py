import math

def evaluate_slater_integral(zeta_i, zeta_j):
    """
    Evaluates the integral <phi_i| 1/r |phi_j> for 1s Slater-type orbitals.

    Args:
        zeta_i (float): The orbital exponent for phi_i.
        zeta_j (float): The orbital exponent for phi_j.
    """
    print(f"Evaluating the integral for zeta_i = {zeta_i} and zeta_j = {zeta_j}")

    if zeta_i <= 0 or zeta_j <= 0:
        print("Orbital exponents (zeta) must be positive.")
        return

    # Check for the special case where i = j
    if zeta_i == zeta_j:
        zeta = zeta_i
        result = zeta
        print("This is a diagonal element (i=j), where zeta_i = zeta_j = zeta.")
        print("The analytical formula simplifies to: <phi|1/r|phi> = zeta")
        print(f"Result = {zeta}")
    else:
        # Use the general formula for off-diagonal elements
        numerator = 4 * (zeta_i * zeta_j)**1.5
        denominator = (zeta_i + zeta_j)**2
        result = numerator / denominator
        
        print("This is an off-diagonal element (i != j).")
        print("The formula is: 4 * (zeta_i * zeta_j)^(3/2) / (zeta_i + zeta_j)^2")
        print(f"<phi_i|1/r|phi_j> = (4 * ({zeta_i} * {zeta_j})^1.5) / ({zeta_i} + {zeta_j})^2")
        print(f"                 = {numerator} / {denominator}")
        print(f"                 = {result}")

    print("-" * 40)
    return result

if __name__ == '__main__':
    # Example 1: Ground state of a Hydrogen atom.
    # For a hydrogenic atom, zeta = Z/n. For H, Z=1, n=1, so zeta=1.
    # This is a diagonal element <phi_1s|1/r|phi_1s>.
    print("Example 1: Hydrogen atom ground state (zeta = 1.0)")
    evaluate_slater_integral(1.0, 1.0)

    # Example 2: A diagonal element for a different exponent,
    # e.g., for a Helium atom using an effective charge (Z_eff ~ 1.69)
    print("Example 2: A different diagonal element (zeta = 1.69)")
    evaluate_slater_integral(1.69, 1.69)
    
    # Example 3: An off-diagonal element between two different 1s STOs.
    print("Example 3: An off-diagonal element (zeta_i=1.0, zeta_j=1.24)")
    evaluate_slater_integral(1.0, 1.24)
