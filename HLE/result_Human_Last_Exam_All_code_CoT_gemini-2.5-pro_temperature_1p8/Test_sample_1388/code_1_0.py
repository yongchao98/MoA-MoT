import sympy as sp

def solve_pde_h_t():
    """
    This function defines the symbolic components of the solution H(t) and prints
    the final expression.
    """
    
    # Define the mathematical symbols using sympy
    # ||u_0||_L2: The L2 norm of the initial condition u_0
    # ||u_0||_L1: The L1 norm of the initial condition u_0, which is the conserved mass M0
    # h(t): The given function h(t) = integral from 0 to t of the L-infinity norm of the spatial derivative of u
    u0_L2_norm_sq = sp.Symbol('||u_0||_L2^2')
    u0_L1_norm = sp.Symbol('||u_0||_L1')
    h_t = sp.Symbol('h(t)')

    # The conserved mass M_0 is the L1 norm of the initial condition u_0
    M0 = u0_L1_norm

    # Numerator of the fraction inside the square root
    # This corresponds to 2 * (1 - e^(-M0)) * h(t) from our derivation
    numerator = 2 * (1 - sp.exp(-M0)) * h_t
    
    # Denominator of the fraction
    denominator = u0_L2_norm_sq

    # The full expression for H(t) is the square root of (1 + numerator/denominator)
    H_t = sp.sqrt(1 + numerator / denominator)

    # Print the derived expression for H(t)
    print("The explicit form of H(t) is:")
    # Using pprint for a more readable mathematical output
    sp.pprint(H_t, use_unicode=True)

solve_pde_h_t()
