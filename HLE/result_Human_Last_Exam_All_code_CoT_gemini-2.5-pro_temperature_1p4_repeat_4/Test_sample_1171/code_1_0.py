def pi_k_S_even_Q_is_zero(k, n_even):
    """
    Checks if pi_k(S^n) tensor Q is zero for an even n.
    For S^{2m}, non-zero rational homotopy groups are at k=2m, 4m-1, etc.
    For our range, only the first two matter.
    """
    m = n_even / 2
    if k == n_even or k == 2 * n_even - 1:
        return False
    return True

def pi_k_CP2_Q_is_zero(k):
    """
    Checks if pi_k(CP^2) tensor Q is zero.
    Non-zero for k=2 and k=5.
    """
    if k == 2 or k == 5:
        return False
    return True

def solve():
    """
    Finds k in {1..9} for which pi_k(S^4 v CP^2) tensor Q is zero.
    This happens if and only if both pi_k(S^4) tensor Q and pi_k(CP^2) tensor Q are zero.
    """
    vanishing_k = []
    n_sphere = 4
    
    for k in range(1, 10):
        # Check if the rational homotopy group of S^4 is zero for this k
        s4_is_zero = pi_k_S_even_Q_is_zero(k, n_sphere)
        
        # Check if the rational homotopy group of CP^2 is zero for this k
        cp2_is_zero = pi_k_CP2_Q_is_zero(k)
        
        if s4_is_zero and cp2_is_zero:
            vanishing_k.append(k)
            
    # Print the result in the required format
    print(",".join(map(str, vanishing_k)))

solve()