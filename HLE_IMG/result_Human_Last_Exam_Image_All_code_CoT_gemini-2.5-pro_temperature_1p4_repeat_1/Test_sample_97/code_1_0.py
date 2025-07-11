import math

def solve():
    # Step 1 & 2: Find relationship between r_w and y_w from bottom row geometry
    # Given yellow circle center Cy = (4, 0.5) and radius r_y = 0.5.
    # It's tangent to two white circles (radius r_w, center y-coord y_w).
    # White circles are tangent to each other. Their centers are Cw1 = (4-r_w, y_w) and Cw2 = (4+r_w, y_w).
    # The distance between Cy and Cw1 must be r_y + r_w = 0.5 + r_w.
    # (r_w + 0.5)^2 = ( (4) - (4-r_w) )^2 + ( (y_w) - (0.5) )^2
    # r_w^2 + r_w + 0.25 = r_w^2 + y_w^2 - y_w + 0.25
    # This simplifies to: r_w = y_w^2 - y_w = y_w * (y_w - 1)
    
    # Step 4: Find r_w and y_w using the 'multiple of 0.5' constraint and x-coordinate data.
    r_w = 2.0 # From the given coordinate x=4, the x of the first yellow circle is 2*r_w. So 2*r_w=4 => r_w=2.
    
    # Now find y_w from r_w = y_w * (y_w-1) => 2 = y_w^2 - y_w => y_w^2 - y_w - 2 = 0
    # Factoring: (y_w - 2)(y_w + 1) = 0. Since y_w > 0, y_w = 2.0
    y_w = 2.0
    r_y = 0.5

    # --- Question 1: Gap between yellow and white circles? ---
    
    # Center of yellow circle
    c_y_x, c_y_y = 4.0, 0.5
    # Center of adjacent white circle
    c_w1_x, c_w1_y = c_y_x - r_w, y_w # (4-2, 2) = (2, 2)
    
    # Distance between centers of yellow and white circles
    dist_yw_sq = (c_y_x - c_w1_x)**2 + (c_y_y - c_w1_y)**2
    dist_yw = math.sqrt(dist_yw_sq)
    
    # Sum of radii
    radii_sum_yw = r_y + r_w
    
    # Compare distance to sum of radii
    gap_yw = dist_yw - radii_sum_yw
    ans1 = 'Y' if gap_yw > 1e-9 else 'N'
    
    print("--- Question 1: Gap between yellow and white circles? ---")
    print(f"Radius of yellow circle (r_y): {r_y}")
    print(f"Radius of white circle (r_w): {r_w}")
    print(f"Center of yellow circle: ({c_y_x}, {c_y_y})")
    print(f"Center of adjacent white circle: ({c_w1_x}, {c_w1_y})")
    print(f"Distance between centers = sqrt(({c_y_x} - {c_w1_x})^2 + ({c_y_y} - {c_w1_y})^2) = {dist_yw}")
    print(f"Sum of radii = {r_y} + {r_w} = {radii_sum_yw}")
    print(f"Result: The distance equals the sum of radii, so there is no gap. Answer: {ans1}")
    print("-" * 20)

    # --- Question 2: Gap between white circles in first and second row? ---
    
    # Step 3: Use the green rectangle to find the vertical separation of white circle rows.
    # "No gap" means the green rectangle is tangent to the circles above and below.
    # This forces the vertical distance between row centerlines to be 2*r_w.
    vertical_dist_centers = 2 * r_w
    horizontal_dist_centers = r_w
    
    # Distance between centers of a bottom-row circle and middle-row circle
    dist_ww_sq = horizontal_dist_centers**2 + vertical_dist_centers**2
    dist_ww = math.sqrt(dist_ww_sq)
    
    # Sum of radii
    radii_sum_ww = r_w + r_w
    
    # Compare
    gap_ww = dist_ww - radii_sum_ww
    ans2 = 'Y' if gap_ww > 1e-9 else 'N'

    print("--- Question 2: Gap between white circles in different rows? ---")
    print(f"Radius of white circle (r_w): {r_w}")
    print(f"Vertical distance between row centerlines (from green rect): 2 * {r_w} = {vertical_dist_centers}")
    print(f"Horizontal distance between centers: {horizontal_dist_centers}")
    print(f"Distance between centers = sqrt({horizontal_dist_centers}^2 + {vertical_dist_centers}^2) = {dist_ww}")
    print(f"Sum of radii = {r_w} + {r_w} = {radii_sum_ww}")
    print(f"Result: The distance is greater than the sum of radii, so there is a gap. Answer: {ans2}")
    print("-" * 20)
    
    print(f"Final Answer: {ans1}{ans2}")
    
solve()