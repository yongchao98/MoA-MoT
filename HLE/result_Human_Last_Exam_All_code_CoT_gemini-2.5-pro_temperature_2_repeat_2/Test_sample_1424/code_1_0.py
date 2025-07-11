import cmath

def calculate_bare_greens_function():
    """
    Calculates and displays the bare Green's function for various
    single-particle energies (epsilon_k).
    """
    # Probe energy (in arbitrary units)
    omega = 5.0
    
    # Small positive infinitesimal for causality
    delta = 0.1
    
    # List of single-particle energy eigenvalues to test
    epsilon_k_values = [4.0, 4.8, 5.0, 5.2, 6.0]
    
    print(f"Calculating the Bare Green's function G_0(k, w) = 1 / (w - e_k + i*d)")
    print(f"Using fixed probe energy w = {omega} and infinitesimal d = {delta}\n")
    
    for e_k in epsilon_k_values:
        # The term in the denominator
        denominator = (omega - e_k) + 1j * delta
        
        # Calculate the Green's function
        # Using 1.0 to ensure floating-point division
        g_0 = 1.0 / denominator
        
        # Print the equation with all numbers and the result
        print(f"For e_k = {e_k}:")
        print(f"  G_0 = 1 / ({omega} - {e_k} + {delta}j)")
        print(f"  G_0 = 1 / {denominator}")
        # Using .2f for cleaner formatting of the complex number result
        print(f"  G_0 = {g_0:.2f}\n")
        
# Execute the function
calculate_bare_greens_function()