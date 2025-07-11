import sys
import math

def solve():
    """
    Solves the chemical bonding simulation problem using dynamic programming with bitmasking.
    """

    def dist_points(p1, p2):
        """Calculates Euclidean distance between two 3D points."""
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

    def dist_point_to_segment(p, a, b):
        """Calculates the shortest distance from a point P to a line segment AB."""
        vec_ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        vec_ap = (p[0] - a[0], p[1] - a[1], p[2] - a[2])
        
        len_sq_ab = vec_ab[0]**2 + vec_ab[1]**2 + vec_ab[2]**2
        
        if len_sq_ab == 0.0:
            return dist_points(p, a)
            
        dot_ap_ab = vec_ap[0] * vec_ab[0] + vec_ap[1] * vec_ab[1] + vec_ap[2] * vec_ab[2]
        t = dot_ap_ab / len_sq_ab
        
        if t < 0.0:
            return dist_points(p, a)
        elif t > 1.0:
            return dist_points(p, b)
        else:
            proj_x = a[0] + t * vec_ab[0]
            proj_y = a[1] + t * vec_ab[1]
            proj_z = a[2] + t * vec_ab[2]
            projection = (proj_x, proj_y, proj_z)
            return dist_points(p, projection)

    lines = sys.stdin.readlines()
    if not lines:
        return
    n = int(lines[0])
    points = [tuple(map(float, line.split())) for line in lines[1:n + 1]]

    if n <= 1:
        print("-1")
        return

    infinity = float('inf')
    
    # --- Step 1: Pre-compute cost for all possible single bonds ---
    memo_cost = [infinity] * (1 << n)
    indices_map = {}
    popcount = [0] * (1 << n)

    for i in range(1, 1 << n):
        popcount[i] = popcount[i >> 1] + (i & 1)
        if popcount[i] >= 2:
            indices = [j for j in range(n) if (i >> j) & 1]
            indices_map[i] = indices

    for mask in range(1, 1 << n):
        if popcount[mask] < 2:
            continue
        
        indices = indices_map[mask]
        min_bond_cost_for_mask = infinity
        
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p_i_idx = indices[i]
                p_j_idx = indices[j]
                
                current_bond_cost = dist_points(points[p_i_idx], points[p_j_idx])
                
                for k in range(len(indices)):
                    if k != i and k != j:
                        p_k_idx = indices[k]
                        current_bond_cost += dist_point_to_segment(points[p_k_idx], points[p_i_idx], points[p_j_idx])
                
                min_bond_cost_for_mask = min(min_bond_cost_for_mask, current_bond_cost)
        memo_cost[mask] = min_bond_cost_for_mask

    # --- Step 2: Dynamic Programming with Path Reconstruction ---
    dp = [infinity] * (1 << n)
    path = [0] * (1 << n)
    dp[0] = 0

    for mask in range(1, 1 << n):
        sub = mask
        while sub > 0:
            group_cost = memo_cost[sub]
            if group_cost != infinity:
                remaining_mask = mask ^ sub
                if dp[remaining_mask] != infinity:
                    new_cost = dp[remaining_mask] + group_cost
                    if new_cost < dp[mask]:
                        dp[mask] = new_cost
                        path[mask] = sub
            sub = (sub - 1) & mask

    final_cost = dp[(1 << n) - 1]

    if final_cost == infinity:
        print("-1")
        return
        
    # --- Step 3: Backtrack and print the equation ---
    costs = []
    current_mask = (1 << n) - 1
    while current_mask > 0:
        group_mask = path[current_mask]
        costs.append(memo_cost[group_mask])
        current_mask ^= group_mask

    equation_parts = [f"{c:.4f}" for c in sorted(costs)]
    print(" + ".join(equation_parts), f"= {final_cost:.4f}")

solve()