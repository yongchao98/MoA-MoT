import sys

def solve():
    """
    Calculates the minimal number of generators for the fundamental group of the
    double branched cover of S^4 over the 5-twist-spun knot of the trefoil.
    """
    
    # Step 1: Define the Alexander polynomial of the trefoil knot.
    # Delta(t) = t - 1 + t^-1
    
    # Step 2: Evaluate the polynomial at t = -1 to find the order of the
    # first homology group of the 3D branched cover (Sigma_2(T)).
    t = -1
    # Using integer division for t^-1 as 1/t
    delta_val = t - 1 + (1 / t)
    
    # The order of the group is the absolute value of this evaluation.
    # The fundamental group of Sigma_2(T) is Z_order.
    order = abs(int(delta_val))
    
    # Step 3: The fundamental group of the 4D manifold is the same, Z_order.
    
    # Step 4: Find the minimal number of generators for Z_order.
    # For any cyclic group Z_n with n > 1, the minimal number of generators is 1.
    if order > 1:
        min_generators = 1
    elif order == 1:
        min_generators = 0 # The trivial group
    else: # order == 0, group is Z
        min_generators = 1
    
    print("The Alexander polynomial of the trefoil knot is Delta(t) = t - 1 + 1/t.")
    print(f"Evaluating at t = -1 gives: Delta(-1) = {t} - 1 + {1/t:.0f} = {int(delta_val)}.")
    print(f"The order of the fundamental group of the associated 3-manifold is |{int(delta_val)}| = {order}.")
    print(f"The fundamental group of the manifold in question is Z_{order}.")
    print(f"The minimal number of generators for Z_{order} is {min_generators}.")

solve()