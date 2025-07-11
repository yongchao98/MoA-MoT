import math

def solve_joule_heat():
    """
    This function explains the derivation and prints the final formula for the Joule heat.
    The problem assumes a sphere's potential is held constant at V as it shrinks from radius a to 0.
    """
    
    # Symbolic representation of variables
    V = "V"
    a = "a"
    k = "k"
    epsilon_0 = "epsilon_0"
    pi = "pi"
    
    print("Based on the principle of energy conservation, the total Joule heat (H) dissipated can be determined.")
    print("The final result is derived from the energy balance of the electrostatic system.")
    print("\nThe formula for the total Joule heat (H) dissipated is:")
    
    # Final formula in terms of k
    final_formula_k = f"H = 2 * {V}^2 * {a} / {k}"
    print(f"\nIn terms of Coulomb's constant, k: \n{final_formula_k}")
    
    # Final formula in terms of epsilon_0
    final_formula_epsilon = f"H = 8 * {pi} * {epsilon_0} * {a} * {V}^2"
    print(f"\nIn terms of the permittivity of free space, epsilon_0:")
    print(final_formula_epsilon)
    
    # Per the instruction "output each number in the final equation"
    print("\nFinal symbolic equation:")
    print("H = 8 * pi * epsilon_0 * a * V^2")

# Execute the function to display the solution
solve_joule_heat()