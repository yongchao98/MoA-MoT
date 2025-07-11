import sympy

def solve_christoffel_symbols():
    """
    This function calculates the Christoffel symbols for the Schwarzschild metric,
    prints each non-zero component's formula, and counts the total number of non-zero components.
    """
    # Step 1: Define the symbolic variables.
    # We use natural units (c=G=1). Coordinates are t, r, theta, phi.
    t, r, theta, phi, M = sympy.symbols('t r theta phi M', real=True)
    coords = [t, r, theta, phi]
    coord_names = ['t', 'r', 'theta', 'phi']

    # Step 2: Define the Schwarzschild metric tensor (g_mu_nu).
    f = 1 - 2 * M / r
    g_cov = sympy.zeros(4)  # Covariant metric tensor g_mu_nu
    g_cov[0, 0] = -f
    g_cov[1, 1] = 1 / f
    g_cov[2, 2] = r**2
    g_cov[3, 3] = r**2 * sympy.sin(theta)**2

    # Step 3: Calculate the inverse metric tensor (g^mu^nu).
    g_con = g_cov.inv()

    # Step 4: Calculate all Christoffel symbols of the second kind.
    # We store the non-zero symbols and their calculated formulas in a list.
    non_zero_symbols = []

    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                # Calculate Gamma^rho_{mu,nu} using the formula.
                sum_val = sympy.sympify(0)
                for sigma in range(4):
                    # Partial derivatives of the metric tensor components
                    term1 = sympy.diff(g_cov[sigma, mu], coords[nu])
                    term2 = sympy.diff(g_cov[sigma, nu], coords[mu])
                    term3 = sympy.diff(g_cov[mu, nu], coords[sigma])
                    
                    sum_val += g_con[rho, sigma] * (term1 + term2 - term3)

                # Final symbol value after simplification.
                christoffel_symbol = sympy.simplify(0.5 * sum_val)

                # If the symbol is not zero, add it to our list for printing.
                if christoffel_symbol != 0:
                    non_zero_symbols.append((rho, mu, nu, christoffel_symbol))
    
    # Step 5: Print the formulas for each non-zero symbol and the total count.
    print("The non-zero Christoffel symbols and their formulas are:")
    for rho, mu, nu, val in non_zero_symbols:
        # We print the symbolic formula for each non-zero component.
        print(f"Gamma^{coord_names[rho]}_({coord_names[mu]},{coord_names[nu]}) = {str(val)}")
        
    print(f"\nTotal number of non-zero Christoffel symbols: {len(non_zero_symbols)}")

if __name__ == '__main__':
    solve_christoffel_symbols()