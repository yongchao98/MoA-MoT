import math

def solve_lindhard_function():
    """
    Determines and explains the numerical value of the Lindhard polarization
    function at zero external frequency and zero momentum transfer.
    """
    
    print("### Solving for the Lindhard Polarization Function Π(k=0, ω=0) ###\n")
    
    print("1. The Physical Relationship:")
    print("The Lindhard polarization function, Π(k, ω), describes the response of an electron gas to a perturbation.")
    print("In the static (ω=0) and long-wavelength (k=0) limit, it is directly related to the")
    print("density of states at the Fermi level, D(E_F), by the fundamental equation:")
    print("  Π(0, 0) = -D(E_F)\n")

    print("2. Dependence on System Parameters:")
    print("For a three-dimensional electron gas, the density of states D(E_F) is not a universal constant.")
    print("It depends on system properties such as the electron density (n) or the Fermi energy (E_F).")
    print("The formula is: D(E_F) = 3*n / (2*E_F).")
    print("Therefore, Π(0, 0) itself does not have a single, universal numerical value.\n")

    print("3. Obtaining a Universal Numerical Value through Normalization:")
    print("The request for a single 'numerical value' implies we should find a dimensionless quantity.")
    print("The natural way to do this is to normalize the Lindhard function by the quantity that sets its scale, which is D(E_F).")
    print("We will evaluate the normalized polarization function: Π(0, 0) / D(E_F).\n")
    
    print("4. Final Calculation:")
    print("Using the relationship from step 1, we can write the equation for the normalized value:")
    
    # In our equation, the result is derived from a ratio of symbolic quantities.
    # The numerator can be represented by -1 (in units of D(E_F))
    # and the denominator by 1 (in units of D(E_F)).
    numerator = -1
    denominator = 1
    
    # The result of the division is the final numerical value.
    final_value = numerator / denominator
    
    # Print the equation representing the calculation.
    # The equation is effectively: (-D(E_F)) / (D(E_F)) = -1
    print(f"  Normalized Π(0, 0) = (-D(E_F)) / (D(E_F))")
    print(f"  The equation simplifies to: {numerator} / {denominator} = {final_value}")

    print("\n----------------------------------------------------------------------")
    print("The numerical value of the Lindhard polarization function at k=0 and ω=0, when normalized by the density of states at the Fermi level, is:")
    print(f"{final_value}")
    print("----------------------------------------------------------------------")


if __name__ == "__main__":
    solve_lindhard_function()
