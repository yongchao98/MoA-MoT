import sys

def solve_bose_light_quanta():
    """
    This function explains and presents the equilibrium values (formulas) for
    mean energy and entropy for a Bose gas of light quanta (photons), as
    derived from the principles of statistical mechanics and large deviation theory.
    """

    print("### Equilibrium Energy and Entropy for Light Quanta (Bose Gas) ###")

    print("\n--- Theoretical Result: The Bose-Einstein Distribution ---")
    print("The starting point is the equilibrium distribution, which results from maximizing entropy.")
    print("The average number of photons, <n_i>, in a state with energy ε_i is:")
    print("\n<n_i> = 1 / (exp(β * ε_i) - 1)\n")
    print("Where:")
    print("  - ε_i: Energy of the i-th state.")
    print("  - β  : Inverse temperature, defined as 1 / (k_B * T).")
    print("  - k_B: Boltzmann's constant.")
    print("  - T  : Absolute temperature.")

    print("\n\n--- 1. Equilibrium Mean Energy (U) ---")
    print("The total mean energy U of the system is the sum of the energies of all photons.")
    print("This is calculated by summing the energy of each state multiplied by its average occupation number.")
    print("\nFinal Equation for Mean Energy:")

    # Printing the formula by breaking it down into its constituent parts
    # as per the request to "output each number in the final equation".
    # Here, '1' is the only explicit number in the formula.
    print("U = ", end="")
    print("Σ_i [ ", end="")
    print("ε_i", end="")
    print(" / ", end="")
    print("(", "exp(β * ε_i)", "-", "1", ")", "]")
    print("\n")
    print("Explanation of components:")
    print("  - U        : Total mean energy.")
    print("  - Σ_i      : Sum over all energy states i.")
    print("  - ε_i      : Energy of the i-th state.")
    print("  - exp(...) : The exponential function.")
    print("  - β        : The inverse temperature, 1/(k_B*T).")
    print("  - 1        : The number one, characteristic of the Bose-Einstein distribution's form.")


    print("\n\n--- 2. Equilibrium Entropy (S) ---")
    print("The entropy S is derived from the thermodynamic relation involving the partition function Z:")
    print("S = k_B * (ln(Z) + β*U). For a Bose gas with zero chemical potential, ln(Z) = -Σ_i ln(1 - exp(-β*ε_i)).")
    print("Substituting the expressions for U and ln(Z) yields the formula for entropy.")
    print("\nFinal Equation for Entropy:")

    # Printing the entropy formula by breaking it down.
    # Here, '1' appears twice.
    print("S = ", end="")
    print("k_B * Σ_i [ ", end="")
    print("(", "(β * ε_i)", "/", "(exp(β * ε_i) - 1)", ")", end="")
    print(" - ", end="")
    print("ln(", "1", "-", "exp(-β * ε_i)", ")", "]")
    print("\n")
    print("Explanation of components:")
    print("  - S        : Total entropy of the system.")
    print("  - k_B      : Boltzmann's constant.")
    print("  - Σ_i      : Sum over all energy states i.")
    print("  - β        : The inverse temperature, 1/(k_B*T).")
    print("  - ε_i      : Energy of the i-th state.")
    print("  - ln(...)  : The natural logarithm.")
    print("  - exp(...) : The exponential function.")
    print("  - 1        : The number one.")

if __name__ == '__main__':
    solve_bose_light_quanta()
