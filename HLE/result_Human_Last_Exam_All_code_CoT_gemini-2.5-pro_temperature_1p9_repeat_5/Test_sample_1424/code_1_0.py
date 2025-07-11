def print_bare_greens_function(epsilon_k, omega):
    """
    Prints the functional form of the bare Green's function G_0
    for a given single-particle energy epsilon_k and frequency omega.

    In the Feynman path integral formalism, the bare Green's function G_0(k, omega)
    is derived from the non-interacting action. It represents the propagator of a
    free particle and is given by the inverse of the operator in the quadratic
    part of the action in frequency-momentum space.

    Args:
        epsilon_k (float): The single-particle energy eigenvalue.
        omega (float): The frequency.
    """
    
    print("The bare Green's function G_0(k, omega) is functionally dependent on the single-particle energy, epsilon_k, as its inverse.")
    print("The general formula for the time-ordered propagator is:")
    print("  G_0(k, omega) = 1 / (omega - epsilon_k + i*delta)")
    print("\nwhere 'delta' is a positive infinitesimal ensuring correct causal structure.")
    print("----------------------------------------------------------------------")
    
    # Illustrate the formula with the provided numerical values
    print(f"For a single-particle energy epsilon_k = {epsilon_k} and frequency omega = {omega}, the expression is:")
    
    # Construct and print the final equation string
    # This shows the structure with the specific numbers.
    print(f"\n  G_0 = 1 / ({omega} - {epsilon_k} + i*delta)")


if __name__ == '__main__':
    # Example values for a single-particle state
    # These could be any physically relevant numbers, e.g., in units of eV.
    example_epsilon_k = 2.5
    example_omega = 3.1
    
    print_bare_greens_function(example_epsilon_k, example_omega)