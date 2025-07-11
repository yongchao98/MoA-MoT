import math
import sys

def solve():
    """
    Solves the chemical bonding simulation problem by finding the minimum cost
    to connect all compounds into bonds.
    """
    # Set a high value for infinity
    INF = float('inf')

    # Helper functions for vector operations
    def v_sub(p2, p1):
        """Subtracts vector p1 from p2."""
        return (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])

    def dot_product(v1, v2):
        """Calculates the dot product of two vectors."""
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]

    def distance(p1, p2):
        """Calculates the Euclidean distance between two points."""
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

    def distance_point_segment(p, a, b):
        """Calculates the shortest distance from a point P to a line segment AB."""
        ab = v_sub(b, a)
        ap = v_sub(p, a)
        
        dot_ab_ab = dot_product(ab, ab)
        if dot_ab_ab == 0:
            return distance(p, a)

        t = dot_product(ap, ab) / dot_ab_ab
        if t < 0.0:
            return distance(p, a)
        elif t > 1.0:
            return distance(p, b)
        else:
            projection = (a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2])
            return distance(p, projection)

    # Hardcoded input data from the test case
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

    if N < 2:
        print(-1 if N == 1 else "0.0000")
        return

    # Step 1: Pre-compute cost[mask] for forming a single bond
    cost = [INF] * (1 << N)
    for mask in range(1, 1 << N):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue

        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                p1, p2 = points[p1_idx], points[p2_idx]
                
                current_bond_cost = distance(p1, p2)
                for k_idx in indices:
                    if k_idx != p1_idx and k_idx != p2_idx:
                        pk = points[k_idx]
                        current_bond_cost += distance_point_segment(pk, p1, p2)
                
                cost[mask] = min(cost[mask], current_bond_cost)

    # Step 2: DP with backtracking to find the optimal partition
    dp = [INF] * (1 << N)
    choice = [-1] * (1 << N)
    dp[0] = 0

    for mask in range(1, 1 << N):
        # Fix the group containing the first point in the mask to ensure each partition is generated once
        first_point_idx = (mask & -mask).bit_length() - 1
        
        # Iterate over all subsets of the rest, to form the group with the first point
        submask_of_rest = mask ^ (1 << first_point_idx)
        s = submask_of_rest
        while True:
            group_mask = s | (1 << first_point_idx)
            rem_mask = mask ^ group_mask

            if cost[group_mask] != INF and dp[rem_mask] != INF:
                current_total_cost = dp[rem_mask] + cost[group_mask]
                if current_total_cost < dp[mask]:
                    dp[mask] = current_total_cost
                    choice[mask] = group_mask

            if s == 0:
                break
            s = (s - 1) & submask_of_rest

    final_cost = dp[(1 << N) - 1]

    if final_cost == INF:
        print(-1)
        return

    # Step 3: Backtrack to find the bonds in the optimal partition and print the "equation"
    bond_costs = []
    current_mask = (1 << N) - 1
    while current_mask > 0:
        best_submask = choice[current_mask]
        bond_costs.append(cost[best_submask])
        current_mask ^= best_submask
        
    bond_costs.sort(reverse=True)
    
    # Step 4: Print the final result per instructions
    formatted_costs = [f"{c:.4f}" for c in bond_costs]
    print(f"The equation for the minimum cost is: {' + '.join(formatted_costs)} = {final_cost:.4f}")
    print("\nFinal Answer (minimum cost):")
    print(f"{final_cost:.4f}")

solve()
<<<492.2030>>>