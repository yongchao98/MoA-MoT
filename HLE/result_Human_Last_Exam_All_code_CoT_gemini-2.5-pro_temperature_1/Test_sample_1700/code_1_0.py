import math
import sys

def solve():
    """
    Solves the chemical bonding simulation problem using dynamic programming with bitmasking.
    """
    # A very large number to represent infinity
    INF = float('inf')

    # Helper class for 3D points and vector operations
    class Point:
        def __init__(self, x, y, z):
            self.x, self.y, self.z = x, y, z

        def __sub__(self, other):
            return Point(self.x - other.x, self.y - other.y, self.z - other.z)

        def dot(self, other):
            return self.x * other.x + self.y * other.y + self.z * other.z

        def magnitude_sq(self):
            return self.dot(self)

        def magnitude(self):
            return math.sqrt(self.magnitude_sq())

    def dist_points(p1, p2):
        """Calculates the Euclidean distance between two points."""
        return (p1 - p2).magnitude()

    def dist_point_to_segment(p, a, b):
        """Calculates the shortest distance from point p to line segment ab."""
        ab = b - a
        ap = p - a
        
        ab_mag_sq = ab.magnitude_sq()
        
        # If the segment is just a point
        if ab_mag_sq == 0:
            return ap.magnitude()
            
        # Project p onto the line defined by a and b to find the parameter t
        t = ab.dot(ap) / ab_mag_sq
        
        if t < 0.0:
            # Closest point on segment is a
            return ap.magnitude()
        elif t > 1.0:
            # Closest point on segment is b
            bp = p - b
            return bp.magnitude()
        else:
            # Projection falls on the segment
            projection_vector = Point(t * ab.x, t * ab.y, t * ab.z)
            projection_point = Point(a.x + projection_vector.x, a.y + projection_vector.y, a.z + projection_vector.z)
            return (p - projection_point).magnitude()

    # --- Main Logic ---
    try:
        lines = sys.stdin.readlines()
        if not lines:
            return
            
        n_str = lines[0].strip()
        if not n_str:
            return
        n = int(n_str)
        
        points = []
        for i in range(1, n + 1):
            x, y, z = map(float, lines[i].split())
            points.append(Point(x, y, z))
    except (IOError, ValueError):
        print("Invalid input format.")
        return

    if n == 1:
        print(-1.0)
        return

    # --- Step 1: Pre-calculate cost for each possible group (subset) ---
    num_masks = 1 << n
    group_cost = [INF] * num_masks

    for mask in range(1, num_masks):
        indices = [i for i in range(n) if (mask >> i) & 1]

        if len(indices) < 2:
            continue

        min_cost_for_this_group = INF
        
        # Iterate through all pairs (i, j) to form the bridge
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                p1, p2 = points[p1_idx], points[p2_idx]
                
                # Bridge cost
                current_cost = dist_points(p1, p2)
                
                # Add cost of connecting other compounds to the bridge
                for k in range(len(indices)):
                    if k != i and k != j:
                        pk_idx = indices[k]
                        pk = points[pk_idx]
                        current_cost += dist_point_to_segment(pk, p1, p2)
                
                min_cost_for_this_group = min(min_cost_for_this_group, current_cost)
        
        group_cost[mask] = min_cost_for_this_group

    # --- Step 2: DP to find minimum cost partition ---
    dp = [INF] * num_masks
    dp[0] = 0.0

    for mask in range(1, num_masks):
        submask = mask
        while submask > 0:
            if dp[mask ^ submask] != INF and group_cost[submask] != INF:
                current_partition_cost = dp[mask ^ submask] + group_cost[submask]
                dp[mask] = min(dp[mask], current_partition_cost)
            submask = (submask - 1) & mask
            
    result = dp[num_masks - 1]
    
    if result == INF:
        print(-1.0)
    else:
        print(f"{result:.4f}")

solve()