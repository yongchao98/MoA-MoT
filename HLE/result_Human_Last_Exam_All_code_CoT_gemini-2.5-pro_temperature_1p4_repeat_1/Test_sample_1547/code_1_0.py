import math

def demonstrate_transience_condition(d, num_terms):
    """
    This function demonstrates that the constructed set A is transient
    by showing that the series for the expected number of visits converges.
    
    The condition for transience is that E_0[N(A)] < infinity.
    E_0[N(A)] = sum_{y in A} G(0, y)
    
    We approximate this sum by sum_k |A_k| * G(0, c_k).
    
    For d-dimensional space:
    - Cube k has center c_k = (3^k, 0, ..., 0).
    - Distance |c_k| = 3^k.
    - Cube k has radius L_k = k.
    - Size of boundary |A_k| is approximately 2*d*(2*L_k+1)^(d-1).
    - Green's function G(0, c_k) is proportional to 1 / |c_k|^(d-2).
    
    So the k-th term of the series is proportional to k^(d-1) / (3^k)^(d-2).
    """
    
    if d < 3:
        print("The argument is for d >= 3.")
        return

    print(f"Demonstrating transience for d={d}.")
    print("We check the convergence of the series whose terms are proportional to |A_k| / |c_k|^(d-2).")
    print("-" * 50)
    
    current_sum = 0.0
    # We start from k=1.
    for k in range(1, num_terms + 1):
        # Cube radius L_k = k
        L_k = k
        # Number of points on the boundary |A_k|
        # For simplicity, we use the approximation 2*d*(2*L_k)^(d-1)
        # The exact number is more complex, but this captures the polynomial growth.
        size_A_k = 2 * d * (2 * L_k)**(d - 1)
        
        # Center of the cube c_k
        dist_c_k = 3**k
        
        # Proportionality constant for Green's function is ignored as we only care about convergence
        green_function_approx = 1 / (dist_c_k)**(d - 2)
        
        term = size_A_k * green_function_approx
        current_sum += term
        
        print(f"k = {k:2d}: Term = {term:10.8f}, Partial Sum = {current_sum:10.8f}")
        
    print("-" * 50)
    print(f"The partial sum over {num_terms} terms is {current_sum:.8f}.")
    print("The terms decrease rapidly, suggesting the series converges.")
    print("This supports the argument that the constructed set A is transient.")

# Let's run the demonstration for d=3 and 15 terms.
demonstrate_transience_condition(d=3, num_terms=15)
