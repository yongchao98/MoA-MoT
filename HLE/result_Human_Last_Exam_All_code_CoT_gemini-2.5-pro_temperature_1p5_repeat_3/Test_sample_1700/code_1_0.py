import math

def solve_chemical_bonding():
    """
    Calculates the minimum cost to connect N compounds using a dynamic programming
    approach with bitmasking.
    """

    # --- Helper functions for geometric calculations ---

    def dist_sq(p1, p2):
        """Calculates the squared Euclidean distance between two 3D points."""
        return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2

    def dist_points(p1, p2):
        """Calculates the Euclidean distance between two 3D points."""
        return math.sqrt(dist_sq(p1, p2))

    def dist_point_to_segment(p, a, b):
        """Calculates the shortest distance from point p to the line segment ab."""
        ab_sq = dist_sq(a, b)
        if ab_sq < 1e-12:  # Treat as the same point if very close
            return dist_points(p, a)

        # Project vector ap onto ab, find parameter t
        dot_product = ((p[0] - a[0]) * (b[0] - a[0]) +
                       (p[1] - a[1]) * (b[1] - a[1]) +
                       (p[2] - a[2]) * (b[2] - a[2]))
        t = dot_product / ab_sq

        if t < 0.0:
            return dist_points(p, a)
        if t > 1.0:
            return dist_points(p, b)
        
        # Projection is on the segment. Use Pythagorean theorem for perpendicular distance.
        ap_sq = dist_sq(p, a)
        dist_squared_to_line = ap_sq - (dot_product * t)
        
        # Clamp to 0 to handle potential floating point inaccuracies
        return math.sqrt(max(0, dist_squared_to_line))

    # --- Input Processing ---

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
    91.725218 58.914653 -26.53959
    39.454154 41.937457 0.591171
    """
    lines = input_data.strip().split('\n')
    points = [tuple(map(float, line.split())) for line in lines]
    
    if N < 2:
        if N <= 0: return
        print(-1)
        return

    # --- Step 1: Pre-compute costs for single bonds ---
    
    cost = [float('inf')] * (1 << N)
    for mask in range(1, 1 << N):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue

        min_bond_cost = float('inf')
        # Iterate through all pairs in the subset to act as the bridge
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                
                bridge_cost = dist_points(points[p1_idx], points[p2_idx])
                attachment_cost = 0.0
                
                # Calculate cost for other points to attach to this bridge
                for k in range(len(indices)):
                    if k != i and k != j:
                        p_other_idx = indices[k]
                        attachment_cost += dist_point_to_segment(
                            points[p_other_idx], points[p1_idx], points[p2_idx]
                        )
                
                total_cost = bridge_cost + attachment_cost
                min_bond_cost = min(min_bond_cost, total_cost)
        
        cost[mask] = min_bond_cost

    # --- Step 2: Dynamic Programming with path tracking ---
    
    dp = [float('inf')] * (1 << N)
    parent = [0] * (1 << N)
    dp[0] = 0.0
    
    for mask in range(1, 1 << N):
        # Option 1: The current set forms a single bond
        dp[mask] = cost[mask]
        parent[mask] = mask  # Sentinel value for a single bond group
        
        # Option 2: The current set is a union of two smaller, disjoint sets
        # Iterate through all proper submasks of the current mask
        submask = (mask - 1) & mask
        while submask > 0:
            other_submask = mask ^ submask
            
            combined_cost = dp[submask] + dp[other_submask]
            if combined_cost < dp[mask]:
                dp[mask] = combined_cost
                parent[mask] = submask  # Store the submask that created this split
                
            submask = (submask - 1) & mask
    
    final_cost = dp[(1 << N) - 1]

    # --- Step 3: Reconstruct solution and Print ---

    if final_cost == float('inf'):
        print(-1)
    else:
        # Reconstruct the partition to find individual bond costs
        bond_costs = []
        q = [(1 << N) - 1]
        while q:
            current_mask = q.pop(0)
            p_mask = parent[current_mask]
            
            if p_mask == current_mask:
                # This is a fundamental bond group
                bond_costs.append(cost[current_mask])
            else:
                # This group was split, so process its children
                q.append(p_mask)
                q.append(current_mask ^ p_mask)
        
        equation_parts = [f"{c:.4f}" for c in sorted(bond_costs, reverse=True)]
        equation_str = " + ".join(equation_parts)
        print(f"{equation_str} = {final_cost:.4f}")

        # Final answer as required by the format
        print(f"<<<{final_cost:.4f}>>>")

solve_chemical_bonding()