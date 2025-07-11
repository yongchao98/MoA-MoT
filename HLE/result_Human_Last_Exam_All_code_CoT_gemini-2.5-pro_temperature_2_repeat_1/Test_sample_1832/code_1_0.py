import sympy as sp

def solve_christoffel_symbols():
    """
    Calculates the number of non-zero Christoffel symbols for the Schwarzschild metric.
    """
    # 1. Define symbolic variables for coordinates and mass M.
    # We set G=c=1 for simplicity.
    t, r, theta, phi, M = sp.symbols('t r theta phi M', real=True, positive=True)

    # 2. Define the coordinates
    coords = [t, r, theta, phi]

    # 3. Define the covariant Schwarzschild metric tensor g_(\mu\nu)
    g_cov = sp.zeros(4)
    g_cov[0, 0] = -(1 - 2*M/r)
    g_cov[1, 1] = (1 - 2*M/r)**(-1)
    g_cov[2, 2] = r**2
    g_cov[3, 3] = r**2 * sp.sin(theta)**2

    # 4. Calculate the contravariant (inverse) metric tensor g^(\mu\nu)
    # For a diagonal metric, this is the element-wise reciprocal.
    g_con = g_cov.applyfunc(lambda x: 1/x if x != 0 else 0)
    
    # 5. Pre-calculate the derivatives of the metric tensor with respect to all coordinates.
    # dg[k,i,j] = \partial_k g_{ij}
    dg = sp.Array([sp.diff(g_cov, c) for c in coords])
    
    # 6. Initialize a list to hold the non-zero symbols and their calculated values
    non_zero_symbols = []

    # 7. Iterate through all 64 possible Christoffel symbols
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                # Calculate the Christoffel symbol Gamma^rho_{mu,nu} using the formula
                # Note that for a diagonal inverse metric, the sum over sigma collapses to a single term sigma=rho.
                term = sp.S.Zero
                # Sum over sigma is implicit in g_con[rho, rho]
                term = g_con[rho, rho] * (dg[nu, rho, mu] + dg[mu, rho, nu] - dg[rho, mu, nu])
                
                # The final Christoffel symbol
                christoffel_symbol = sp.simplify(sp.S(1)/2 * term)

                # 8. Check if the simplified symbol is non-zero
                if christoffel_symbol != sp.S.Zero:
                    non_zero_symbols.append(christoffel_symbol)

    # 9. Print the final results as requested
    count = len(non_zero_symbols)
    print(f"The total number of non-zero Christoffel symbols is: {count}")
    
    # "Remember in the final code you still need to output each number in the final equation!"
    # This prints an equation like 1 + 1 + ... + 1 = total count.
    count_str = " + ".join(["1" for _ in range(count)])
    if not count_str: # handles case of zero
        count_str = "0"
        
    print(f"The calculation is: {count_str} = {count}")

solve_christoffel_symbols()