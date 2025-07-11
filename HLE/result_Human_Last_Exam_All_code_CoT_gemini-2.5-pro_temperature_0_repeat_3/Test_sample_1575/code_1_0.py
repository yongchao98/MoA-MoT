import math

def solve_reversal_moves():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given adjacent and non-adjacent swap rules.
    """
    N = 100
    MOD = 5

    # 1. There are N/MOD = 20 elements in each of the MOD=5 classes.
    elements_per_class = N // MOD

    # 2. Determine the required migration of elements between classes.
    # An element at initial position p_i = k-1 (1-indexed k) moves to p_f = N-k.
    # Initial class C_i = (k-1) % MOD.
    # Final class C_f = (N-k) % MOD.
    # If C_i = C, then k = q*MOD + C + 1.
    # C_f = (N - (q*MOD + C + 1)) % MOD = (-C - 1) % MOD.
    # This gives the mapping from a starting class to a destination class.
    
    # 3. Calculate the flows between classes based on shortest paths on a 5-node cycle.
    # flows_forward[i] stores the number of elements crossing from class i to (i+1).
    # flows_backward[i] stores the number of elements crossing from class (i+1) to i.
    flows_forward = [0] * MOD
    flows_backward = [0] * MOD

    for start_class in range(MOD):
        end_class = (-start_class - 1 + MOD) % MOD
        
        # Elements starting and ending in the same class require no inter-class moves.
        if start_class == end_class:
            continue
        
        count = elements_per_class
        
        # Determine the shortest path on the cycle graph (forward or backward).
        dist_forward = (end_class - start_class + MOD) % MOD
        dist_backward = (start_class - end_class + MOD) % MOD

        if dist_forward <= dist_backward:
            # Move forward: start -> start+1 -> ... -> end
            curr = start_class
            for _ in range(dist_forward):
                flows_forward[curr] += count
                curr = (curr + 1) % MOD
        else:
            # Move backward: start -> start-1 -> ... -> end
            curr = start_class
            for _ in range(dist_backward):
                prev = (curr - 1 + MOD) % MOD
                # This flow is from curr to prev, which is a backward flow over the (prev, curr) boundary.
                flows_backward[prev] += count
                curr = prev

    # 4. Calculate the number of swaps required at each boundary.
    # A swap at boundary (i, i+1) facilitates one forward and one backward move.
    # So, the number of swaps is the maximum of the flows in either direction.
    swaps_per_boundary = [0] * MOD
    for i in range(MOD):
        swaps_per_boundary[i] = max(flows_forward[i], flows_backward[i])

    # 5. The total number of moves is the sum of swaps across all boundaries.
    total_moves = sum(swaps_per_boundary)

    # Print the final calculation
    equation_parts = [str(s) for s in swaps_per_boundary]
    equation = " + ".join(equation_parts)
    print("The minimum number of moves is the sum of swaps required at each class boundary.")
    print(f"Calculation: {equation} = {total_moves}")

solve_reversal_moves()