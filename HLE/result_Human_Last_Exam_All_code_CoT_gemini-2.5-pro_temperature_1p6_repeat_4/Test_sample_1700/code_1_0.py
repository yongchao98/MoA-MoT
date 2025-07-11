import math
import sys

def solve_chemical_bonding():
    """
    This script calculates the minimum cost to connect N compounds by partitioning
    them into bonds, using dynamic programming over subsets.
    """
    try:
        lines = sys.stdin.readlines()
        if not lines:
            return
        n_str = lines[0].strip()
        if not n_str:
            return
        n = int(n_str)
        
        if n == 0:
            print("0.0000")
            print("<<<0.0000>>>")
            return
        
        # Per problem constraints, a bond requires at least 2 compounds.
        # It's impossible to connect just one compound, so we print -1.
        if n == 1:
            print("-1")
            print("<<<-1>>>")
            return
            
        points_str = lines[1:n+1]
        points = []
        for line in points_str:
            x, y, z = map(float, line.strip().split())
            points.append((x, y, z))

    except (IOError, ValueError, IndexError):
        print("Error: Invalid input format.", file=sys.stderr)
        return

    # --- Step 1: Geometric Helper Functions ---
    def dist_sq(p1, p2):
        """Calculates the squared Euclidean distance between two points."""
        return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2

    def dist_point_to_segment(p, a, b):
        """Calculates the shortest distance from a point p to the line segment ab."""
        ab_sq = dist_sq(a, b)
        if ab_sq == 0:  # a and b are the same point
            return math.sqrt(dist_sq(p, a))
        
        dot_product = ((p[0] - a[0]) * (b[0] - a[0]) +
                       (p[1] - a[1]) * (b[1] - a[1]) +
                       (p[2] - a[2]) * (b[2] - a[2]))
        
        t = max(0.0, min(1.0, dot_product / ab_sq))
        
        proj = (a[0] + t * (b[0] - a[0]),
                a[1] + t * (b[1] - a[1]),
                a[2] + t * (b[2] - a[2]))
        
        return math.sqrt(dist_sq(p, proj))

    # --- Step 2: Pre-compute single bond costs for all subsets ---
    bond_cost = [float('inf')] * (1 << n)

    for mask in range(1, 1 << n):
        indices = [i for i in range(n) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue

        min_current_bond_cost = float('inf')
        
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                p1, p2 = points[p1_idx], points[p2_idx]

                current_cost = math.sqrt(dist_sq(p1, p2))
                
                for k in range(len(indices)):
                    if k != i and k != j:
                        pk_idx = indices[k]
                        pk = points[pk_idx]
                        current_cost += dist_point_to_segment(pk, p1, p2)
                
                min_current_bond_cost = min(min_current_bond_cost, current_cost)
        
        bond_cost[mask] = min_current_bond_cost

    # --- Step 3: Dynamic Programming for Optimal Set Partition ---
    dp = [float('inf')] * (1 << n)
    dp_path = [-1] * (1 << n) # To store submask for backtracking
    dp[0] = 0.0

    for mask in range(1, 1 << n):
        sub = mask
        while sub > 0:
            cost_of_sub = bond_cost[sub]
            cost_of_rest = dp[mask ^ sub]
            
            if cost_of_sub != float('inf') and cost_of_rest != float('inf'):
                new_total_cost = cost_of_rest + cost_of_sub
                if new_total_cost < dp[mask]:
                    dp[mask] = new_total_cost
                    dp_path[mask] = sub # Store this submask as part of the optimal path
            
            sub = (sub - 1) & mask
    
    final_cost = dp[(1 << n) - 1]

    # --- Step 4: Output the results ---
    if math.isinf(final_cost):
        print("-1")
        print("<<<-1>>>")
    else:
        # Backtrack to find the bonds that form the optimal solution
        optimal_bonds_masks = []
        current_mask = (1 << n) - 1
        while current_mask > 0:
            best_sub = dp_path[current_mask]
            optimal_bonds_masks.append(best_sub)
            current_mask ^= best_sub

        bond_costs = [bond_cost[mask] for mask in optimal_bonds_masks]
        
        # Sort costs for a consistent output format
        equation_parts = [f"{c:.4f}" for c in sorted(bond_costs, reverse=True)]
        equation_str = " + ".join(equation_parts)

        print(f"{equation_str} = {final_cost:.4f}")
        print(f"<<<{final_cost:.4f}>>>")

if __name__ == "__main__":
    solve_chemical_bonding()