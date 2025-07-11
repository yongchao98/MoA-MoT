import sys
import math

class Point:
    """A helper class for 3D points and vector operations."""
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z

    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y, self.z - other.z)

    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y, self.z + other.z)

    def __rmul__(self, scalar):
        return Point(self.x * scalar, self.y * scalar, self.z * scalar)

    def dot(self, other):
        return self.x * other.x + self.y * other.y + self.z * other.z

    def length_sq(self):
        return self.dot(self)

    def length(self):
        return math.sqrt(self.length_sq())

def distance_3d(p1, p2):
    """Calculates the Euclidean distance between two 3D points."""
    return (p1 - p2).length()

def dist_point_to_segment_3d(p, a, b):
    """Calculates the shortest distance from a point p to a line segment ab in 3D."""
    ab = b - a
    ap = p - a
    
    len_sq_ab = ab.length_sq()
    if len_sq_ab < 1e-12:  # a and b are essentially the same point
        return ap.length()

    # Project vector ap onto ab, find parameter t
    t = ap.dot(ab) / len_sq_ab
    
    if t < 0.0:
        # Closest point on the line is outside the segment, near a
        return ap.length()
    elif t > 1.0:
        # Closest point on the line is outside the segment, near b
        bp = p - b
        return bp.length()
    else:
        # Projection is on the segment
        projection = a + t * ab
        return (p - projection).length()

def solve_simulation():
    """
    Reads compound data, calculates, and prints the minimum bonding cost.
    """
    try:
        lines = sys.stdin.read().strip().split('\n')
        if not lines or not lines[0]:
            # This handles cases where input is piped but empty
            print("Invalid input: Empty file.")
            return

        n = int(lines[0])
        points = []
        if len(lines) < n + 1:
            print(f"Invalid input: Expected {n} coordinate lines, but got {len(lines)-1}.")
            return
            
        for i in range(1, n + 1):
            points.append(Point(*map(float, lines[i].split())))

    except (IOError, ValueError, IndexError) as e:
        print(f"Error processing input: {e}. Please ensure input is formatted correctly.")
        return

    # A bond requires at least 2 compounds.
    if n < 2:
        print("-1")
        return

    num_masks = 1 << n
    
    # Step 1: Pre-calculate the minimum cost for every subset to form a single bond
    bond_cost = [float('inf')] * num_masks
    
    for mask in range(1, num_masks):
        indices = [i for i in range(n) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue
            
        # Iterate over all pairs in the subset to act as the bridge
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                
                # Cost for this bridge
                current_cost = distance_3d(points[p1_idx], points[p2_idx])
                
                # Add costs for other compounds in the subset connecting to the bridge
                for k in range(len(indices)):
                    if k != i and k != j:
                        pk_idx = indices[k]
                        current_cost += dist_point_to_segment_3d(points[pk_idx], points[p1_idx], points[p2_idx])
                
                bond_cost[mask] = min(bond_cost[mask], current_cost)

    # Step 2: Use DP with bitmasking to find the minimum cost partitioning
    dp = [float('inf')] * num_masks
    path = [-1] * num_masks  # To reconstruct the final partition for the equation
    dp[0] = 0
    
    for mask in range(1, num_masks):
        submask = mask
        while submask > 0:
            remaining_mask = mask ^ submask
            cost = bond_cost[submask] + dp[remaining_mask]
            
            if cost < dp[mask]:
                dp[mask] = cost
                path[mask] = submask  # Store which submask gave the optimal cost
                
            submask = (submask - 1) & mask

    final_cost = dp[num_masks - 1]
    
    if final_cost == float('inf'):
        print("-1")
        return

    # Step 3: Reconstruct the partition and print the result as an equation
    partition_masks = []
    current_mask = num_masks - 1
    while current_mask > 0:
        best_submask = path[current_mask]
        partition_masks.append(best_submask)
        current_mask ^= best_submask

    partition_costs = [bond_cost[m] for m in partition_masks]
    equation_parts = [f"{c:.4f}" for c in sorted(partition_costs)]
    
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {final_cost:.4f}")

if __name__ == '__main__':
    solve_simulation()