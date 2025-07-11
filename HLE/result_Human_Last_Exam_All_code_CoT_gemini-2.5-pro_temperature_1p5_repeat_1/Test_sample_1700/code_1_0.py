import math
import sys

def solve():
    """
    Solves the chemical bonding problem using dynamic programming with bitmasking.
    """
    
    # --- Helper functions for 3D geometry ---
    
    def dist(p1, p2):
        """Calculates Euclidean distance between two 3D points."""
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

    def dist_point_segment(p, a, b):
        """Calculates the shortest distance from point p to the line segment ab."""
        ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        ap = (p[0] - a[0], p[1] - a[1], p[2] - a[2])

        ab_len_sq = ab[0]**2 + ab[1]**2 + ab[2]**2
        
        if ab_len_sq == 0.0:
            return dist(p, a)
            
        dot_product = ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]
        t = dot_product / ab_len_sq

        if t < 0.0:
            return dist(p, a)
        elif t > 1.0:
            return dist(p, b)
        else:
            proj_point = (a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2])
            return dist(p, proj_point)

    # --- Main logic ---

    try:
        lines = sys.stdin.read().strip().split('\n')
        if not lines or not lines[0]:
            return
            
        n = int(lines[0])
        points = [tuple(map(float, line.split())) for line in lines[1:n+1]]

    except (IOError, ValueError):
        print("Invalid input format.", file=sys.stderr)
        return
        
    if n < 2:
        print(-1)
        return

    num_subsets = 1 << n
    
    # 1. Pre-compute cost of forming a single bond for every subset
    bond_cost = [float('inf')] * num_subsets

    for mask in range(1, num_subsets):
        indices = [i for i in range(n) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue

        min_s_cost = float('inf')
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                p1, p2 = points[p1_idx], points[p2_idx]
                
                current_cost = dist(p1, p2) # Cost of the bridge
                
                for k in range(len(indices)):
                    if k != i and k != j:
                        p3_idx = indices[k]
                        p3 = points[p3_idx]
                        current_cost += dist_point_segment(p3, p1, p2)
                
                min_s_cost = min(min_s_cost, current_cost)
        
        bond_cost[mask] = min_s_cost

    # 2. DP to find the minimum cost partition
    dp = [float('inf')] * num_subsets
    path = [-1] * num_subsets  # To reconstruct the partition
    dp[0] = 0.0

    for mask in range(1, num_subsets):
        # Iterate over all non-empty submasks of the current mask
        submask = mask
        while submask > 0:
            remaining_mask = mask ^ submask
            cost_of_this_partition = bond_cost[submask] + dp[remaining_mask]
            
            if cost_of_this_partition < dp[mask]:
                dp[mask] = cost_of_this_partition
                path[mask] = submask # Store the submask that forms one bond
                
            submask = (submask - 1) & mask

    # 3. Reconstruct path and print the final equation
    total_cost = dp[num_subsets - 1]
    
    if total_cost == float('inf'):
         print(-1)
         return
         
    costs = []
    current_mask = num_subsets - 1
    while current_mask > 0:
        bond_mask = path[current_mask]
        costs.append(bond_cost[bond_mask])
        current_mask ^= bond_mask

    formatted_costs = [f"{c:.4f}" for c in costs]
    total_cost_formatted = f"{total_cost:.4f}"
    
    print(f"{' + '.join(formatted_costs)} = {total_cost_formatted}")


solve()