def check_potts_submodularity():
    """
    Checks the submodularity condition for the Potts model energy on a single edge.

    The Gibbs measure is log-supermodular if the energy H is submodular, i.e.,
    H(xi) + H(eta) >= H(xi v eta) + H(xi ^ eta).
    For Potts model H is a sum over edges of H_e, so we check one edge.
    H_e(sigma) = -I(sigma(x) == sigma(y)).
    The condition becomes:
    -I(xi(x)==xi(y)) - I(eta(x)==eta(y)) >= -I(max...) - I(min...)
    or equivalently:
    I(xi(x)==xi(y)) + I(eta(x)==eta(y)) <= I(max...) + I(min...)

    We show a counterexample for q=3.
    """
    print("Checking the submodularity condition for the Potts energy function.")
    print("The condition for a single edge is: I(xi_x=xi_y) + I(eta_x=eta_y) <= I(max(xi_x,eta_x)=max(xi_y,eta_y)) + I(min(xi_x,eta_x)=min(xi_y,eta_y))\n")

    # Counterexample configuration for q=3
    xi = (2, 2)
    eta = (1, 3)

    print(f"Let's choose q=3 and states from {{1, 2, 3}}.")
    print(f"Consider two configurations on a single edge (x,y):")
    print(f"  xi = (xi(x), xi(y)) = {xi}")
    print(f"  eta = (eta(x), eta(y)) = {eta}\n")

    # Compute coordinate-wise max (join) and min (meet)
    xi_v_eta = (max(xi[0], eta[0]), max(xi[1], eta[1]))
    xi_m_eta = (min(xi[0], eta[0]), min(xi[1], eta[1]))

    print("The derived configurations are:")
    print(f"  xi v eta = {xi_v_eta}")
    print(f"  xi ^ eta = {xi_m_eta}\n")

    # Calculate indicator function values
    I_xi = 1 if xi[0] == xi[1] else 0
    I_eta = 1 if eta[0] == eta[1] else 0
    I_xi_v_eta = 1 if xi_v_eta[0] == xi_v_eta[1] else 0
    I_xi_m_eta = 1 if xi_m_eta[0] == xi_m_eta[1] else 0
    
    # Calculate both sides of the inequality
    lhs = I_xi + I_eta
    rhs = I_xi_v_eta + I_xi_m_eta

    print("Now, we evaluate the inequality:")
    print(f"I(xi(x)=xi(y)) + I(eta(x)=eta(y))  =  {I_xi} + {I_eta}  =  {lhs}")
    print(f"I(join) + I(meet)             =  {I_xi_v_eta} + {I_xi_m_eta}  =  {rhs}")
    
    print(f"\nThe required inequality is {lhs} <= {rhs}.")
    
    if lhs <= rhs:
        print("This is true. The condition holds for this example.")
    else:
        print("This is FALSE.")
        print("\nBecause the submodularity condition fails, the standard proof of the positive correlation property does not apply.")
        print("This failure on a single edge (a graph with max degree d=1) indicates that d>=1 is not a correct answer,")
        print("as the property must hold for all q>=2, beta>=0, and all graphs with max degree d.")
        print("The largest d for which the property is guaranteed is d=0.")


check_potts_submodularity()