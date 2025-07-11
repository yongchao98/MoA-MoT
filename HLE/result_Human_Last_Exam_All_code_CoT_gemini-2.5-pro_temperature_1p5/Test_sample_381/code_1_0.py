import math

def calculate_upper_bound(N, M, lamb, deltas):
    """
    Calculates the upper bound for ||B * Q_{0,M}||_infinity.

    Args:
        N (int): The number of nodes in the graph.
        M (int): The number of layers/time steps in the product.
        lamb (float): The uniform spectral norm bound for projected stochastic matrices (must be < 1).
        deltas (list[float]): A list of M+1 delta values, where delta_k = ||D^(k) - I_N||_infinity.
    """
    if len(deltas) != M + 1:
        raise ValueError(f"The length of deltas must be M+1 = {M+1}.")

    # Start with the lambda^(M+1) term, corresponding to a_{-1}=1
    # a_M <= lambda^(M+1) * a_{-1} + sum_{k=0 to M} lambda^(M-k) * delta_k
    a_M_bound = lamb ** (M + 1)
    
    # Calculate the summation part
    summation_term = 0
    for k in range(M + 1):
        term = (lamb ** (M - k)) * deltas[k]
        summation_term += term
    
    a_M_bound += summation_term
    
    # The final bound is sqrt(N) * a_M_bound
    final_bound = math.sqrt(N) * a_M_bound
    
    # Print the equation with numerical values
    print("Upper bound calculation:")
    print(f"||B Q_{{0,M}}||_inf <= sqrt(N) * (lambda^(M+1) + sum_{{k=0 to M}} lambda^(M-k) * delta_k)")
    
    # Print sqrt(N) part
    print(f"                 = sqrt({N}) * (", end="")
    
    # Print lambda^(M+1) part
    lambda_term_val = lamb ** (M + 1)
    print(f"{lamb}^({M}+1)", end="")

    # Print summation part
    full_sum_val = 0
    for k in range(M + 1):
        term_val = (lamb ** (M - k)) * deltas[k]
        full_sum_val += term_val
        print(f" + {lamb}^({M-k})*{deltas[k]}", end="")
    print(")")
    
    # Print intermediate computed values
    print(f"                 = {math.sqrt(N):.4f} * ({lambda_term_val:.4f} + {full_sum_val:.4f})")
    
    # Print the final result
    print(f"                 = {final_bound:.4f}")

    return final_bound

# Example usage with some chosen parameters
N_nodes = 50
M_layers = 5
lambda_bound = 0.95
# delta_k values for k = 0, 1, 2, 3, 4, 5
delta_values = [0.1, 0.08, 0.06, 0.04, 0.02, 0.01]

calculate_upper_bound(N_nodes, M_layers, lambda_bound, delta_values)