import sympy

def calculate_christoffel_symbols():
    """
    This function calculates the non-zero Christoffel symbols for the Schwarzschild metric,
    prints a breakdown of the components, and outputs the total count.
    """
    # 1. Define symbolic variables
    t, r, theta, phi = sympy.symbols('t r theta phi')
    M, G, c = sympy.symbols('M G c', positive=True)
    
    # Coordinates and their names for printing
    coords = [t, r, theta, phi]
    coord_names = ['t', 'r', 'theta', 'phi']

    # Schwarzschild radius shorthand
    r_s = 2 * G * M / c**2

    # 2. Define the Schwarzschild metric tensor (covariant components)
    g_cov = sympy.Matrix([
        [-(1 - r_s/r) * c**2, 0, 0, 0],
        [0, 1 / (1 - r_s/r), 0, 0],
        [0, 0, r**2, 0],
        [0, 0, 0, r**2 * sympy.sin(theta)**2]
    ])

    # 3. Calculate the inverse metric (contravariant components)
    g_con = g_cov.inv()

    # 4. Calculate derivatives of the metric with respect to all coordinates
    # dg[i, j, k] = d(g_cov[i,j]) / d(coords[k])
    dg = sympy.derive_by_array(g_cov, coords)

    # 5. Calculate and store unique non-zero Christoffel symbols
    # Use a dictionary to store unique symbols, using canonical indices (rho, min(mu,nu), max(mu,nu))
    # to handle symmetry Gamma^rho_munu = Gamma^rho_numu
    unique_symbols = {}
    
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                term = sympy.S.Zero
                for sigma in range(4):
                    term += g_con[rho, sigma] * (dg[sigma, mu, nu] + dg[sigma, nu, mu] - dg[mu, nu, sigma])
                
                christoffel_symbol = sympy.simplify(term / 2)

                if christoffel_symbol != 0:
                    # Canonical key for the unique symbol
                    canonical_key = (rho, min(mu, nu), max(mu, nu))
                    if canonical_key not in unique_symbols:
                         unique_symbols[canonical_key] = {
                             "value": christoffel_symbol,
                             "indices": (rho, mu, nu)
                         }

    # 6. Print the breakdown and count the total components
    print("The non-zero Christoffel symbols for the Schwarzschild metric are calculated as follows:\n")

    total_count = 0
    count_list = []
    
    # Sort by upper index, then lower indices for a structured output
    sorted_keys = sorted(unique_symbols.keys())

    for key in sorted_keys:
        symbol_info = unique_symbols[key]
        rho, mu, nu = symbol_info["indices"]
        
        # Determine number of components for this unique symbol
        num_components = 1 if mu == nu else 2
        
        # Build the name string
        name = f"Gamma^({coord_names[rho]})_({coord_names[mu]},{coord_names[nu]})"
        if mu != nu:
            name += f" = Gamma^({coord_names[rho]})_({coord_names[nu]},{coord_names[mu]})"
        
        print(f"- {name}")
        # print(f"  Value: {symbol_info['value']}") # This would print the full expression
        print(f"  This corresponds to {num_components} non-zero component(s).")
        
        total_count += num_components
        count_list.append(str(num_components))

    print("\n" + "="*50)
    print("Final Calculation:")
    print("The total number of non-zero Christoffel symbols is the sum of these components.")
    print(f"{' + '.join(count_list)} = {total_count}")
    print("="*50)


if __name__ == '__main__':
    calculate_christoffel_symbols()