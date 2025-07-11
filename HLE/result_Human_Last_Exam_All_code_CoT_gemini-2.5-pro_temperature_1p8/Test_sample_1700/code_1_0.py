import sys
import math

class Point:
    """A helper class to represent a 3D point."""
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

def dist_points(p1, p2):
    """Calculates the Euclidean distance between two 3D points."""
    return math.sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2)

def dist_point_to_segment(p, a, b):
    """Calculates the shortest distance from point p to the line segment ab."""
    ab = Point(b.x - a.x, b.y - a.y, b.z - a.z)
    ap = Point(p.x - a.x, p.y - a.y, p.z - a.z)
    
    len_sq_ab = ab.x**2 + ab.y**2 + ab.z**2
    if len_sq_ab == 0:
        return dist_points(p, a)

    dot_product = ap.x * ab.x + ap.y * ab.y + ap.z * ab.z
    t = dot_product / len_sq_ab
    
    if t < 0:
        return dist_points(p, a)
    elif t > 1:
        return dist_points(p, b)
    else:
        projection = Point(a.x + t * ab.x, a.y + t * ab.y, a.z + t * ab.z)
        return dist_points(p, projection)

def solve_reaction():
    """Main function to solve the problem."""
    try:
        lines = sys.stdin.readlines()
        if not lines: return
        n = int(lines[0])
        
        if n < 2:
            print("-1.0000")
            return
            
        points = []
        for i in range(1, n + 1):
            coords = lines[i].split()
            points.append(Point(coords[0], coords[1], coords[2]))
    except (ValueError, IndexError):
        print("Invalid input format.")
        return

    num_masks = 1 << n
    cost = [float('inf')] * num_masks
    
    # Pre-compute cost of forming a single bond for each subset
    for mask in range(num_masks):
        indices = [i for i in range(n) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue

        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                p1, p2 = points[p1_idx], points[p2_idx]
                
                current_bond_cost = dist_points(p1, p2)
                
                for k in range(len(indices)):
                    if k == i or k == j: continue
                    pk_idx = indices[k]
                    pk = points[pk_idx]
                    current_bond_cost += dist_point_to_segment(pk, p1, p2)
                
                cost[mask] = min(cost[mask], current_bond_cost)

    # DP with bitmasking and parent tracking for equation reconstruction
    dp = [float('inf')] * num_masks
    parent = [-1] * num_masks
    dp[0] = 0

    for mask in range(1, num_masks):
        submask = mask
        while submask > 0:
            if cost[submask] != float('inf') and dp[mask ^ submask] != float('inf'):
                new_cost = dp[mask ^ submask] + cost[submask]
                if new_cost < dp[mask]:
                    dp[mask] = new_cost
                    parent[mask] = submask
            submask = (submask - 1) & mask

    final_total_cost = dp[num_masks - 1]

    if final_total_cost == float('inf'):
         print("-1.0000")
         return

    # Reconstruct the partition to print the equation
    partition_costs = []
    current_mask = num_masks - 1
    while current_mask > 0:
        sub_mask_for_partition = parent[current_mask]
        partition_costs.append(cost[sub_mask_for_partition])
        current_mask ^= sub_mask_for_partition
    
    # Sort for a consistent output format
    partition_costs.sort(reverse=True)
    output_parts = [f"{c:.4f}" for c in partition_costs]
    
    equation_str = " + ".join(output_parts)
    print(f"{equation_str} = {final_total_cost:.4f}")

solve_reaction()