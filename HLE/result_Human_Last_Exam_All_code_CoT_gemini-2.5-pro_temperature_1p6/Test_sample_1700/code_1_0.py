import math

class Point:
    """A simple class to represent a point in 3D space."""
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

def dist(p1, p2):
    """Calculates the Euclidean distance between two points."""
    return math.sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2)

def point_segment_dist(p, a, b):
    """
    Calculates the shortest distance from point p to the line segment ab.
    This is done by projecting p onto the infinite line defined by a and b.
    If the projection falls within the segment, the distance is the perpendicular
    distance. Otherwise, the distance is to the closest endpoint (a or b).
    """
    # Squared length of the segment vector.
    ab_dist_sq = (a.x - b.x)**2 + (a.y - b.y)**2 + (a.z - b.z)**2
    if ab_dist_sq == 0:
        return dist(p, a)

    # Project p onto the line using the dot product.
    # t = dot(p - a, b - a) / |b - a|^2
    dot_product = (p.x - a.x) * (b.x - a.x) + \
                  (p.y - a.y) * (b.y - a.y) + \
                  (p.z - a.z) * (b.z - a.z)
    t = dot_product / ab_dist_sq

    if t < 0:
        # Projection is outside the segment, beyond a. Closest point is a.
        return dist(p, a)
    elif t > 1:
        # Projection is outside the segment, beyond b. Closest point is b.
        return dist(p, b)
    else:
        # Projection is on the segment. Calculate distance to the projection point.
        proj_point = Point(a.x + t * (b.x - a.x),
                           a.y + t * (b.y - a.y),
                           a.z + t * (b.z - a.z))
        return dist(p, proj_point)

def solve_chemical_bonding():
    """
    Solves the problem by pre-computing bond costs and using dynamic programming
    on subsets to find the minimum cost partition.
    """
    # Step 1: Parse the input
    N = 13
    input_data = """
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
    points = [Point(*line.split()) for line in input_data.strip().split('\n')]
    
    if N <= 1:
        print(-1)
        return
        
    # Step 2: Pre-compute minimum cost for every possible single bond (subset)
    min_bond_cost = [float('inf')] * (1 << N)
    for mask in range(1, 1 << N):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue
            
        current_mask_min_cost = float('inf')
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                
                bridge_cost = dist(points[p1_idx], points[p2_idx])
                
                attachment_cost = 0.0
                for k in range(len(indices)):
                    if k != i and k != j:
                        pk_idx = indices[k]
                        attachment_cost += point_segment_dist(points[pk_idx], points[p1_idx], points[p2_idx])
                
                current_mask_min_cost = min(current_mask_min_cost, bridge_cost + attachment_cost)
        
        min_bond_cost[mask] = current_mask_min_cost

    # Step 3: Dynamic Programming on subsets
    dp = [float('inf')] * (1 << N)
    dp[0] = 0.0
    
    for mask in range(1, 1 << N):
        submask = mask
        while submask > 0:
            if dp[mask ^ submask] != float('inf') and min_bond_cost[submask] != float('inf'):
                dp[mask] = min(dp[mask], dp[mask ^ submask] + min_bond_cost[submask])
            submask = (submask - 1) & mask

    # Step 4: Output the result
    result = dp[(1 << N) - 1]
    print(f"{result:.4f}")

solve_chemical_bonding()