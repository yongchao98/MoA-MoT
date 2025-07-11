import math

def solve_chemical_reaction():
    """
    Solves the chemical reaction simulation problem using dynamic programming with bitmasking.
    """

    # Using a simple class for Point/Vector operations to keep the code clean.
    class Point:
        def __init__(self, x, y, z):
            self.x, self.y, self.z = x, y, z

        def __sub__(self, other):
            return Point(self.x - other.x, self.y - other.y, self.z - other.z)

        def __add__(self, other):
            return Point(self.x + other.x, self.y + other.y, self.z + other.z)
        
        def __mul__(self, scalar):
            return Point(self.x * scalar, self.y * scalar, self.z * scalar)

        def dot(self, other):
            return self.x * other.x + self.y * other.y + self.z * other.z

        def length_sq(self):
            return self.dot(self)

        def length(self):
            return math.sqrt(self.length_sq())

    # Input data from the problem description
    N = 13
    coords_str = """
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
    
    points = []
    for line in coords_str.strip().split('\n'):
        x, y, z = map(float, line.split())
        points.append(Point(x, y, z))

    if N <= 1:
        print(-1)
        return

    # Helper function for distance between two points
    def dist(p1, p2):
        return (p1 - p2).length()

    # Helper function for shortest distance from a point to a line segment
    def dist_point_to_segment(pk, pi, pj):
        v = pj - pi
        w = pk - pi
        l2 = v.length_sq()
        if l2 == 0.0:  # pi and pj are the same point
            return dist(pk, pi)
        
        # Project pk onto the line containing the segment pi-pj.
        # The parameter t is ((pk - pi) dot (pj - pi)) / |pj - pi|^2
        t = w.dot(v) / l2
        
        # Clamp t to the range [0, 1] to stay on the segment
        t = max(0, min(1, t))
        
        # The closest point on the segment is pi + t * v
        projection = pi + v * t
        return dist(pk, projection)

    # --- Part 1: Pre-compute the cost of a single bond for each subset ---
    cost = [float('inf')] * (1 << N)
    
    for mask in range(1, 1 << N):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue
        
        min_mask_cost = float('inf')
        
        # Iterate through all pairs in the subset to form the bridge
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                idx_i, idx_j = indices[i], indices[j]
                p_i, p_j = points[idx_i], points[idx_j]
                
                # Cost of the bridge itself
                current_bond_cost = dist(p_i, p_j)
                
                # Add cost of connecting other compounds to the bridge
                for k in range(len(indices)):
                    if k == i or k == j:
                        continue
                    idx_k = indices[k]
                    p_k = points[idx_k]
                    current_bond_cost += dist_point_to_segment(p_k, p_i, p_j)
                
                min_mask_cost = min(min_mask_cost, current_bond_cost)
        
        cost[mask] = min_mask_cost

    # --- Part 2: DP to find the minimum cost partition ---
    dp = [float('inf')] * (1 << N)
    dp[0] = 0
    
    for mask in range(1, 1 << N):
        # Option 1: The entire 'mask' set forms a single bond.
        dp[mask] = cost[mask]
        
        # Option 2: The 'mask' set is partitioned into two smaller groups.
        # We iterate through all non-empty proper submasks of 'mask'.
        submask = (mask - 1) & mask
        while submask > 0:
            other_submask = mask ^ submask
            partition_cost = dp[submask] + dp[other_submask]
            if partition_cost < dp[mask]:
                dp[mask] = partition_cost
            submask = (submask - 1) & mask

    # The final answer is the minimum cost for the set of all N points.
    final_cost = dp[(1 << N) - 1]
    
    print(f"{final_cost:.4f}")

solve_chemical_reaction()