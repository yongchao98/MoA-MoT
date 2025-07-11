import math
import sys

class Point:
    """A simple class for 3D points to improve readability."""
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

def dist_sq(p1, p2):
    """Calculates the squared Euclidean distance between two points."""
    return (p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2

def dist(p1, p2):
    """Calculates the Euclidean distance between two points."""
    return math.sqrt(dist_sq(p1, p2))

def shortest_dist_to_segment(p, a, b):
    """
    Calculates the shortest distance from a point 'p' to a line segment 'ab'.
    """
    ab_sq = dist_sq(a, b)
    # If the segment is just a point, return distance to that point.
    if ab_sq == 0:
        return dist(p, a)

    # Project vector (p-a) onto vector (b-a).
    # The parameter t represents the position of the projection on the infinite line.
    dot_product = (p.x - a.x) * (b.x - a.x) + \
                  (p.y - a.y) * (b.y - a.y) + \
                  (p.z - a.z) * (b.z - a.z)
    
    t = dot_product / ab_sq

    # If t < 0, the projection is outside the segment, closest to a.
    if t < 0:
        return dist(p, a)
    # If t > 1, the projection is outside the segment, closest to b.
    elif t > 1:
        return dist(p, b)
    # Otherwise, the projection is on the segment. The distance is the perpendicular distance.
    else:
        # Using cross product for perpendicular distance: |(p-a) x (b-a)| / |b-a|
        # We work with squared magnitudes to avoid an intermediate sqrt.
        pa_x, pa_y, pa_z = p.x - a.x, p.y - a.y, p.z - a.z
        ba_x, ba_y, ba_z = b.x - a.x, b.y - a.y, b.z - a.z
        
        cross_x = pa_y * ba_z - pa_z * ba_y
        cross_y = pa_z * ba_x - pa_x * ba_z
        cross_z = pa_x * ba_y - pa_y * ba_x
        
        cross_mag_sq = cross_x**2 + cross_y**2 + cross_z**2
        
        return math.sqrt(cross_mag_sq / ab_sq)

def solve_chemical_bonds():
    """
    Main function to solve the problem by reading input, calculating costs, and printing the result.
    """
    # Step 1: Read input from a multiline string for this specific problem instance.
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
    try:
        N = int(lines[0])
        points = []
        for i in range(1, N + 1):
            x, y, z = map(float, lines[i].split())
            points.append(Point(x, y, z))
    except (IOError, ValueError):
        # Handle potential empty input or format errors
        N = 0
        points = []

    # Step 2: Handle edge case where N < 2 (no bonds possible)
    if N <= 1:
        print(-1.0)
        return

    # Step 3: Pre-compute the minimum cost for forming a single bond for every subset of points.
    bond_cost = {}  # Using a dictionary: mask -> cost
    for mask in range(1, 1 << N):
        # A bond requires at least 2 points.
        if bin(mask).count('1') < 2:
            continue

        indices = [i for i in range(N) if (mask >> i) & 1]
        num_points_in_mask = len(indices)
        
        min_cost_for_mask = float('inf')
        
        # Iterate through all pairs (i, j) in the subset to form the bridge.
        for i in range(num_points_in_mask):
            for j in range(i + 1, num_points_in_mask):
                p1_idx, p2_idx = indices[i], indices[j]
                
                # Cost of the bridge itself.
                current_cost = dist(points[p1_idx], points[p2_idx])
                
                # Add cost of connecting other points to this bridge.
                for k in range(num_points_in_mask):
                    if k != i and k != j:
                        pk_idx = indices[k]
                        current_cost += shortest_dist_to_segment(points[pk_idx], points[p1_idx], points[p2_idx])
                
                min_cost_for_mask = min(min_cost_for_mask, current_cost)
        
        bond_cost[mask] = min_cost_for_mask

    # Step 4: Use Dynamic Programming to find the minimum cost partition.
    dp = [float('inf')] * (1 << N)
    dp[0] = 0

    for mask in range(1, 1 << N):
        # Iterate through all submasks of the current mask to find the optimal partition.
        submask = mask
        while submask > 0:
            if submask in bond_cost:
                remaining_mask = mask ^ submask
                if dp[remaining_mask] != float('inf'):
                    new_cost = dp[remaining_mask] + bond_cost[submask]
                    dp[mask] = min(dp[mask], new_cost)
            
            submask = (submask - 1) & mask

    # Step 5: The final answer is the DP value for the mask representing all points.
    final_cost = dp[(1 << N) - 1]
    
    # Step 6: Print the final calculated cost.
    print(f"{final_cost:.4f}")

# Execute the solution
solve_chemical_bonds()