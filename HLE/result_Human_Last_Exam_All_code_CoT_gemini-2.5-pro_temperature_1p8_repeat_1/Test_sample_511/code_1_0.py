import math

def find_optimal_box():
    """
    Finds the integer dimensions (l, w, h) of a rectangular box that minimize
    surface area while being able to theoretically contain 27 balls of radius 2.
    """
    # Number of balls and their radius
    num_balls = 27
    radius = 2.0
    
    # Calculate the total volume of the balls
    volume_of_balls = num_balls * (4/3) * math.pi * (radius**3)
    
    # Use the densest packing factor (FCC/HCP packing) to find the minimum required box volume
    # Packing density = pi / (3 * sqrt(2))
    packing_density = math.pi / (3 * math.sqrt(2))
    min_required_volume = volume_of_balls / packing_density

    # Initial box surface area to beat
    initial_sa = 6 * (12**2)
    
    # Initialize variables to store the best found solution
    min_sa = initial_sa
    best_dims = (12, 12, 12)
    
    # Search for optimal integer dimensions l, w, h.
    # We can cap the search range since very long thin boxes will have large surface areas.
    # A limit of 40 for each dimension is more than sufficient.
    search_limit = 40
    # To avoid permutations, we search with l <= w <= h
    for l in range(1, search_limit + 1):
        for w in range(l, search_limit + 1):
            for h in range(w, search_limit + 1):
                volume = l * w * h
                
                # Check if the box has enough volume to theoretically hold the balls
                if volume >= min_required_volume:
                    sa = 2 * (l*w + l*h + w*h)
                    
                    # If this is a better solution, update our minimum
                    if sa < min_sa:
                        min_sa = sa
                        best_dims = (l, w, h)
                        
    # Check if a better box was found
    if min_sa < initial_sa:
        l, w, h = best_dims
        # We need to output the integers representing length, width, height, and surface area
        # The print format is crucial for the final output.
        print(f"{l}:{w}:{h}:{int(min_sa)}")
    else:
        print(0)

find_optimal_box()
