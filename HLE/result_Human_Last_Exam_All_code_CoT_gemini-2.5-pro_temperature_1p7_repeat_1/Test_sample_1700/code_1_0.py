import math

def solve_chemical_bonding():
    """
    Calculates the minimum cost to connect a set of N compounds based on a
    partitioning scheme. The problem is solved using dynamic programming on subsets.
    """
    
    # Step 1: Input Data
    # The problem provides a specific test case with N=13.
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
    points = [tuple(map(float, line.split())) for line in points_str.strip().split('\n')]

    if N <= 1:
        print("-1.0000")
        return

    # Step 2: Geometry Helper Functions
    def dist(p1, p2):
        """Calculates Euclidean distance between two 3D points."""
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

    def dist_point_to_segment(p, a, b):
        """Calculates the shortest distance from a point p to a line segment ab."""
        ab_vec = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        ap_vec = (p[0] - a[0], p[1] - a[1], p[2] - a[2])
        
        mag_sq_ab = ab_vec[0]**2 + ab_vec[1]**2 + ab_vec[2]**2
        if mag_sq_ab == 0:
            return dist(p, a)
        
        dot_product = ap_vec[0] * ab_vec[0] + ap_vec[1] * ab_vec[1] + ap_vec[2] * ab_vec[2]
        t = dot_product / mag_sq_ab

        if t < 0:
            closest_point = a
        elif t > 1:
            closest_point = b
        else:
            closest_point = (a[0] + t * ab_vec[0], a[1] + t * ab_vec[1], a[2] + t * ab_vec[2])
        
        return dist(p, closest_point)

    # Step 3: Pre-compute bond_cost[mask] for all subsets
    bond_cost = [float('inf')] * (1 << N)
    
    for mask in range(1, 1 << N):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue

        min_mask_cost = float('inf')
        num_indices = len(indices)
        for i in range(num_indices):
            for j in range(i + 1, num_indices):
                idx1, idx2 = indices[i], indices[j]
                p1, p2 = points[idx1], points[idx2]
                
                # Cost is bridge distance + sum of distances of other points to the bridge
                current_cost = dist(p1, p2)
                for k in range(num_indices):
                    if k != i and k != j:
                        idx3 = indices[k]
                        p3 = points[idx3]
                        current_cost += dist_point_to_segment(p3, p1, p2)
                
                min_mask_cost = min(min_mask_cost, current_cost)
        bond_cost[mask] = min_mask_cost

    # Step 4: Dynamic Programming to find minimum cost partition
    dp = [float('inf')] * (1 << N)
    predecessor = [-1] * (1 << N) # For backtracking the solution
    dp[0] = 0.0

    for mask in range(1, 1 << N):
        lsb_idx = (mask & -mask).bit_length() - 1
        sub_mask_base = mask ^ (1 << lsb_idx)
        
        # Iterate over all subsets of sub_mask_base to form the other part of the "last group"
        sub = sub_mask_base
        while True:
            # The "last group" `s` is composed of the lsb element and `sub`
            s = sub | (1 << lsb_idx)
            remainder = mask ^ s
            
            if dp[remainder] != float('inf') and bond_cost[s] != float('inf'):
                new_cost = dp[remainder] + bond_cost[s]
                if new_cost < dp[mask]:
                    dp[mask] = new_cost
                    predecessor[mask] = s
            
            if sub == 0:
                break
            sub = (sub - 1) & sub_mask_base

    # Step 5 & 6: Backtrack, Format, and Print the Final Result
    final_cost = dp[(1 << N) - 1]

    if final_cost == float('inf'):
        print("-1.0000")
        return
        
    costs = []
    mask = (1 << N) - 1
    while mask > 0:
        s = predecessor[mask]
        costs.append(bond_cost[s])
        mask = mask ^ s
    
    costs.reverse()
    cost_strs = [f"{c:.4f}" for c in costs]
    equation = " + ".join(cost_strs)
    
    if len(costs) > 1:
        print(f"{equation} = {final_cost:.4f}")
    else:
        print(f"{final_cost:.4f}")

solve_chemical_bonding()