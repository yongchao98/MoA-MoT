import cmath

def calculate_bare_greens_function(omega, epsilon_k, delta=0.1):
    """
    Calculates the bare Green's function G_0(k, omega).

    Args:
        omega (float): The energy/frequency variable.
        epsilon_k (float): The single-particle energy eigenvalue.
        delta (float): A small positive infinitesimal for causality.
    """
    # Calculate the Green's function
    g0 = 1 / ((omega - epsilon_k) + 1j * delta)

    # Print the context and parameters
    print("The formula for the bare Green's function is: G_0(k, ω) = 1 / (ω - ε_k + iδ)")
    print("\nFor the given values:")
    print(f"ω = {omega}")
    print(f"ε_k = {epsilon_k}")
    print(f"δ = {delta}")
    
    # Print the final equation with the numbers substituted in
    print("\nResulting equation:")
    # The following line explicitly prints each number in the equation
    print(f"G_0 = 1 / ({omega} - {epsilon_k} + {delta}j)")
    
    # Print the final numerical result, formatted to 4 decimal places
    print(f"\nCalculated G_0 = {g0.real:.4f}{g0.imag:+.4f}j")

# --- Example Usage ---
# We demonstrate the calculation for a set of parameters.
# Here, we choose a probing energy 'omega' that is very close to the
# particle's energy 'epsilon_k', which should result in a large response (resonance).
print("--- Calculating for a near-resonance case ---")
calculate_bare_greens_function(omega=3.05, epsilon_k=3.0, delta=0.1)
