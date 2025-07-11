import math

def find_optimal_box():
    """
    Analyzes the initial pioneer box and searches for a more efficient one.
    """
    
    # --- Baseline Parameters ---
    ball_diameter = 4
    initial_dim = 12
    
    # Calculate balls and surface area for the initial box
    balls_per_dim_initial = initial_dim // ball_diameter
    min_balls_to_contain = balls_per_dim_initial ** 3
    initial_surface_area = 6 * (initial_dim ** 2)

    # Initialize the search with the current best (the initial box)
    best_config = {
        "l": initial_dim,
        "w": initial_dim,
        "h": initial_dim,
        "sa": initial_surface_area
    }
    
    # Search for a better configuration.
    # The search space for ball counts (nx, ny, nz) doesn't need to be huge.
    # A highly elongated box will have a very large surface area.
    # We search up to the number of balls to contain, which is a reasonable limit.
    search_limit = min_balls_to_contain + 1 
    
    found_better_box = False

    # Iterate through combinations of ball counts (nx, ny, nz)
    for nx in range(1, search_limit):
      # Start subsequent loops from the current value to avoid redundant permutations (e.g., 1x2x3 vs 3x2x1)
      for ny in range(nx, search_limit):
        for nz in range(ny, search_limit):
          
          # Check if this configuration can hold enough balls
          if nx * ny * nz >= min_balls_to_contain:
            
            # For this packing, find the minimal integer box dimensions
            l = nx * ball_diameter
            w = ny * ball_diameter
            h = nz * ball_diameter
            
            # Calculate the surface area
            current_sa = 2 * (l*w + l*h + w*h)
            
            # A more efficient box must have a STRICTLY smaller surface area
            if current_sa < best_config["sa"]:
              found_better_box = True
              # If we find a better one, we would update best_config
              # but our analysis shows this will not happen. For completeness:
              best_config["l"] = l
              best_config["w"] = w
              best_config["h"] = h
              best_config["sa"] = current_sa

    # Based on our analysis, the original 3x3x3 configuration in a 12x12x12 box
    # is optimal. The search will not find a configuration with a smaller surface area.
    # Therefore, the answer is 0.
    
    if found_better_box:
      # This part of the code will not be reached, but is included for logical completeness.
      dims = sorted([best_config['l'], best_config['w'], best_config['h']])
      # Format the output as L:W:H:SurfaceArea
      result = f"{dims[0]}:{dims[1]}:{dims[2]}:{best_config['sa']}"
      print(result)
    else:
      # Print 0 if no better box is found.
      print("0")

find_optimal_box()