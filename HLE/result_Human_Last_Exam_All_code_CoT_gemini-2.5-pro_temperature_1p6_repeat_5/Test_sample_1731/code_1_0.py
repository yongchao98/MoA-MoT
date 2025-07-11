import math

def solve_equilibrium_values():
    """
    Derives and prints the equilibrium values for a photon gas based on
    the principle of maximum entropy from statistical mechanics, which is
    rigorously described by large deviation theorems.
    """
    print("### Equilibrium Values for Bose-Einstein Gas (Photons) ###")
    print("\n--- 1. Principle: Maximum Entropy ---")
    print("The equilibrium state of a physical system is the one that maximizes its entropy (S) subject to constraints.")
    print("This concept, central to statistical mechanics, is formally grounded in large deviation theory.")
    print("For bosons like photons, we maximize the Bose-Einstein entropy formula subject to a constant average energy <U>.")
    print("This is done using the method of Lagrange multipliers, which introduces a parameter β = 1 / (k_B * T).")

    print("\n--- 2. Equilibrium Occupation Number <n_i> ---")
    print("The result of the maximization gives the equilibrium (or average) number of photons, <n_i>,")
    print("for an energy level ε_i with g_i possible states (degeneracy).")
    print("This is the Bose-Einstein distribution for photons (where chemical potential is zero):")

    print("\n<n_i> = g_i / (exp(β * ε_i) - 1)\n")

    print("Components of the equation:")
    print(f"  g_i: The degeneracy, or the number of states at energy level i.")
    print(f"  exp: The exponential function, e.")
    print(f"  β: The inverse temperature parameter, 1 / (k_B * T).")
    print(f"  ε_i: The energy of level i.")
    print(f"  1: The number 1 is subtracted in the denominator, a defining feature of Bose-Einstein statistics.")

    print("\n--- 3. Equilibrium Mean Energy <U> ---")
    print("The total mean energy <U> of the system is the sum of the energies of all photons.")
    print("This is found by summing the energy of each level weighted by its equilibrium occupation number:")

    print("\n<U> = Σ_i [ <n_i> * ε_i ]\n")

    print("Substituting the formula for <n_i>, we get the full expression for mean energy:")

    print("\n<U> = Σ_i [ (g_i * ε_i) / (exp(β * ε_i) - 1) ]\n")
    print("Each term in this sum is composed of:")
    print(f"  Numerator  : g_i * ε_i")
    print(f"  Denominator: exp(β * ε_i) - 1")

    print("\n--- 4. Equilibrium Entropy S ---")
    print("The entropy S is found by substituting the equilibrium <n_i> back into the entropy formula.")
    print("A convenient expression for the entropy is derived from the grand canonical partition function (Ξ):")

    print("\nS = (<U> / T) + k_B * ln(Ξ)\n")

    print("This can be written in terms of β and ε_i as:")

    print("\nS = k_B * [ β * <U> - Σ_i g_i * ln(1 - exp(-β * ε_i)) ]\n")

    print("Components of the final equation:")
    print(f"  k_B: The Boltzmann constant.")
    print(f"  The first part, 'β * <U>', represents the energy component of entropy.")
    print(f"    - Its full form is: β * Σ_i [ (g_i * ε_i) / (exp(β * ε_i) - 1) ]")
    print(f"  The second part is the sum over all levels i of: -g_i * ln(1 - exp(-β * ε_i))")
    print(f"    - The number '1' inside the natural logarithm (ln) originates from the geometric series summation")
    print(f"      used to calculate the partition function for bosons.")

if __name__ == '__main__':
    solve_equilibrium_values()
