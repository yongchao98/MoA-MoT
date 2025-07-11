import math

def calculate_intersection_probability():
    """
    This function calculates the probability that a circle of radius R=6, thrown
    at random onto a grid of mesh 1, intersects exactly 47 cells.
    It uses a Monte Carlo simulation for an approximate solution.
    """
    R = 6.0
    R_sq = R * R
    TARGET_CELLS = 47
    
    def count_intersected_cells(u, v):
        """
        Calculates the number of grid cells intersected by the circle's circumference.
        Args:
            u, v: The coordinates of the circle's center.
        Returns:
            The number of intersected cells.
        """
        # Define a bounding box for the cell indices (k, l) to check.
        k_min = math.floor(u - R)
        k_max = math.ceil(u + R)
        l_min = math.floor(v - R)
        l_max = math.ceil(v + R)
        
        intersected_count = 0
        for k in range(k_min, k_max):
            for l in range(l_min, l_max):
                # The cell is the square region [k, k+1] x [l, l+1].
                
                # 1. Find the squared distance from the center (u,v) to the closest point in the cell.
                closest_x = max(k, min(u, k + 1))
                closest_y = max(l, min(v, l + 1))
                d_min_sq = (u - closest_x)**2 + (v - closest_y)**2
                
                # If the circle is too far to touch the cell, it doesn't intersect.
                if d_min_sq > R_sq:
                    continue
                    
                # 2. Find the squared distance from the center (u,v) to the farthest corner of the cell.
                dx_far = max(u - k, k + 1 - u)
                dy_far = max(v - l, l + 1 - v)
                d_max_sq = dx_far**2 + dy_far**2
    
                # If the entire cell is contained within the circle, the circumference does not intersect it.
                if d_max_sq < R_sq:
                    continue
                
                # If neither of the above, the circumference must intersect the cell.
                intersected_count += 1
                
        return intersected_count

    # Perform a Monte Carlo simulation.
    # N_samples determines the grid size (N_samples x N_samples) for sampling.
    # A larger number gives a better approximation but takes longer to compute.
    N_samples = 400
    
    count_target = 0
    total_samples = N_samples * N_samples

    # We sample centers (u, v) from the square [0, 0.5] x [0, 0.5] due to symmetry.
    for i in range(N_samples):
        u = (i + 0.5) / N_samples * 0.5
        for j in range(N_samples):
            v = (j + 0.5) / N_samples * 0.5
            
            num_cells = count_intersected_cells(u, v)
            
            if num_cells == TARGET_CELLS:
                count_target += 1
                
    probability = count_target / total_samples
    
    print("The simulation was run with the following parameters:")
    print(f"Circle Radius: {R}")
    print(f"Target Number of Intersected Cells: {TARGET_CELLS}")
    print(f"Number of sample points in the simulation: {total_samples}")
    print("\nResult:")
    print(f"The number of times the circle intersected exactly {TARGET_CELLS} cells was: {count_target}")
    print(f"The total number of random throws simulated was: {total_samples}")
    print(f"The estimated probability is the ratio of these two numbers:")
    print(f"P(N=47) = {count_target} / {total_samples} = {probability:.4f}")

# Running the simulation to find the answer.
# Please note that this calculation might take a few moments to complete.
calculate_intersection_probability()