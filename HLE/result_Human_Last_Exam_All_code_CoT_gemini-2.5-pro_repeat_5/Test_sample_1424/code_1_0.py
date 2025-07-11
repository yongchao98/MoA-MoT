import cmath

def calculate_bare_greens_function(omega_n, epsilon_k, mu):
    """
    Calculates the bare Green's function G_0(k, i*omega_n).
    
    Args:
        omega_n (float): The Matsubara frequency.
        epsilon_k (float): The single-particle energy eigenvalue.
        mu (float): The chemical potential.
    
    Returns:
        complex: The value of the Green's function.
    """
    # The term i*omega_n in the denominator
    i_omega_n = 1j * omega_n
    
    # The term for the energy relative to the chemical potential
    xi_k = epsilon_k - mu
    
    # G_0 = 1 / (i*omega_n - xi_k)
    g0 = 1 / (i_omega_n - xi_k)
    
    return g0

def main():
    """
    Demonstrates the functional dependence of G_0 on epsilon_k.
    """
    # We choose the lowest positive fermionic Matsubara frequency for demonstration.
    # omega_n = (2n+1)*pi*T. Let's set T=1/(pi) and n=0, so omega_n = 1.0.
    omega_n = 1.0
    
    # We set the chemical potential to zero for simplicity.
    mu = 0.0
    
    # A range of single-particle energy eigenvalues to test
    epsilon_k_values = [-2.0, -1.0, 0.0, 1.0, 2.0]
    
    print(f"Demonstrating G_0 dependence on epsilon_k for i*omega_n = {1j*omega_n} and mu = {mu}\n")
    print("Formula: G_0(k, i*omega_n) = 1 / (i*omega_n - (epsilon_k - mu))\n")

    for epsilon_k in epsilon_k_values:
        # Calculate the Green's function
        g0_value = calculate_bare_greens_function(omega_n, epsilon_k, mu)
        
        # Print the full equation with the numbers substituted in
        i_omega_n_term = 1j * omega_n
        
        # We format the numbers to make the output clear
        print(f"For epsilon_k = {epsilon_k:5.2f}: G_0 = 1 / ({i_omega_n_term:.2f} - ({epsilon_k:5.2f} - {mu:.2f})) = {g0_value.real:.3f} + {g0_value.imag:.3f}j")

if __name__ == "__main__":
    main()