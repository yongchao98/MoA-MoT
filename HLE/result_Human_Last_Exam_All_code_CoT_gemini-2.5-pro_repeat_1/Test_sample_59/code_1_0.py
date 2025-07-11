import numpy as np

def calculate_sbm_edge_probability():
    """
    Calculates the probability of drawing an edge in a specific
    Stochastic Block Model (SBM), which is a jointly exchangeable graph model.

    The probability 'p' is given by the formula: p = pi^T * P * pi
    where:
    - pi is the vector of community probabilities.
    - P is the matrix of inter- and intra-community edge probabilities.
    """
    # Define the parameters for a 3-community Stochastic Block Model
    
    # pi: Probabilities for a node to be in community 1, 2, or 3.
    pi = np.array([0.5, 0.3, 0.2])

    # P: Symmetric matrix of connection probabilities.
    # P[k,l] is the probability of an edge between a node in community k
    # and a node in community l.
    P = np.array([[0.7, 0.1, 0.05],
                  [0.1, 0.5, 0.2],
                  [0.05, 0.2, 0.4]])

    print("Calculating the edge probability for a Stochastic Block Model with:")
    print(f"Community probabilities (pi) = {pi}")
    print("Connection probability matrix (P):")
    print(P)
    print("\n")
    
    # The formula is p = sum over all k,l of (pi_k * pi_l * P_kl)
    # This can be broken down into diagonal and off-diagonal terms
    # p = sum_k(pi_k^2 * P_kk) + sum_{k<l}(2 * pi_k * pi_l * P_kl)

    equation_parts = []
    
    # Add diagonal terms to the equation string
    for i in range(len(pi)):
        part = f"({pi[i]}**2 * {P[i,i]})"
        equation_parts.append(part)
        
    # Add off-diagonal terms to the equation string
    for i in range(len(pi)):
        for j in range(i + 1, len(pi)):
            part = f"(2 * {pi[i]} * {pi[j]} * {P[i,j]})"
            equation_parts.append(part)

    # Calculate the final numerical result
    prob = pi.T @ P @ pi
    
    # Assemble and print the final equation
    equation_str = "p = " + " + ".join(equation_parts)
    print("The probability equation is:")
    print(equation_str)
    
    # Print the final calculated probability
    print(f"\nThe calculated probability p = {prob:.3f}")

calculate_sbm_edge_probability()