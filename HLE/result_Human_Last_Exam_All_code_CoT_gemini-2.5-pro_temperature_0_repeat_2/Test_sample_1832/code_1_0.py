import sympy

def solve_christoffel_symbols():
    """
    Calculates and counts the non-zero Christoffel symbols for the Schwarzschild metric.
    """
    # Define coordinates and constants (using natural units G=c=1)
    # coords = (t, r, theta, phi)
    t, r, theta, phi = sympy.symbols('t r theta phi', real=True)
    M = sympy.Symbol('M', real=True, positive=True)
    coords = [t, r, theta, phi]

    # Define the Schwarzschild metric tensor g_munu
    # Using r_s = 2M for Schwarzschild radius for simplicity in expressions
    rs = 2 * M
    g = sympy.diag(
        -(1 - rs/r),
        1 / (1 - rs/r),
        r**2,
        r**2 * sympy.sin(theta)**2
    )

    # Calculate the inverse metric tensor g^munu
    g_inv = g.inv()
    # Simplify the inverse metric components
    for i in range(4):
        for j in range(4):
            g_inv[i, j] = sympy.simplify(g_inv[i, j])

    # Calculate the derivatives of the metric tensor with respect to each coordinate
    # dg[k, i, j] will store the partial derivative of g_ij with respect to x^k
    dg = sympy.Array([sympy.diff(g, c) for c in coords])

    # Store non-zero Christoffel symbols
    non_zero_symbols = []
    count = 0

    print("Calculating non-zero Christoffel symbols Γ^ρ_{μν}...")
    # Calculate Christoffel symbols Gamma^rho_{mu,nu}
    # Using the formula: Γ^ρ_{μν} = 1/2 * g^{ρσ} * (∂_μ g_{νσ} + ∂_ν g_{μσ} - ∂_σ g_{μν})
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                gamma = sympy.sympify(0)
                for sigma in range(4):
                    # Accessing derivatives: dg[k, i, j] is ∂_k g_ij
                    term = dg[mu, nu, sigma] + dg[nu, mu, sigma] - dg[sigma, mu, nu]
                    gamma += (1/2) * g_inv[rho, sigma] * term
                
                # Simplify the final expression for the symbol
                gamma_simplified = sympy.simplify(gamma)

                if gamma_simplified != 0:
                    count += 1
                    # Using unicode for greek letters for better readability
                    symbol_str = f"Γ^{rho}_{mu},{nu}"
                    non_zero_symbols.append((symbol_str, gamma_simplified))

    # Print the results
    print("\nThe non-zero Christoffel symbols and their values are:")
    for symbol, value in non_zero_symbols:
        print(f"{symbol} = {value}")
    
    print("\nTo find the total number of non-zero Christoffel symbols, we count each distinct non-zero component Γ^ρ_{μν}.")
    print(f"The total number of non-zero Christoffel symbols is: {count}")

if __name__ == '__main__':
    solve_christoffel_symbols()