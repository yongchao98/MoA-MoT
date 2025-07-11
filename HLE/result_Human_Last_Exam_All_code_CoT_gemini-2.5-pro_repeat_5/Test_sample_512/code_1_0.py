import math

def calculate_container_properties(centers, ball_radius_cm):
    """Calculates the bounding box and surface area for a set of ball centers."""
    if not centers:
        return 0, 0, 0, 0

    min_x = min(c[0] for c in centers)
    max_x = max(c[0] for c in centers)
    min_y = min(c[1] for c in centers)
    max_y = max(c[1] for c in centers)
    min_z = min(c[2] for c in centers)
    max_z = max(c[2] for c in centers)

    # Container dimensions are the span of centers plus a radius on each side.
    L = (max_x - min_x) + 2 * ball_radius_cm
    W = (max_y - min_y) + 2 * ball_radius_cm
    H = (max_z - min_z) + 2 * ball_radius_cm

    surface_area = 2 * (L * W + L * H + W * H)
    return L, W, H, surface_area

def run_analysis():
    """
    Analyzes packing efficiency to find a better container.
    """
    initial_sa = 864.0
    ball_radius = 2.0
    
    print("Initial Container: 12.0x12.0x12.0 box")
    print(f"Initial ball capacity: 27")
    print(f"Initial surface area: {initial_sa} cm^2")
    print("-" * 30)
    print("Searching for a more efficient packing strategy...")
    print("Strategy: A discretized BCC-like packing.")
    print("In-plane spacing: 3.0 cm, Layer separation: 3.5 cm.")
    print("-" * 30)

    # This packing alternates between 3x3 layers and 2x2 layers.
    # Base coordinates for the grids (we will offset them)
    grid_3x3 = [(x*3.0, y*3.0) for x in range(3) for y in range(3)] # 9 balls
    grid_2x2 = [(1.5 + x*3.0, 1.5 + y*3.0) for x in range(2) for y in range(2)] # 4 balls
    
    centers = []
    current_z = ball_radius # Start z-position for the first layer's centers
    
    packing_results = []

    # Pack layer by layer up to 5 layers
    for i in range(5):
        if i % 2 == 0: # Add a 3x3 layer
            layer_centers = [(x + ball_radius, y + ball_radius, current_z) for x, y in grid_3x3]
        else: # Add a 2x2 layer
            layer_centers = [(x + ball_radius, y + ball_radius, current_z) for x, y in grid_2x2]
        
        centers.extend(layer_centers)
        
        L, W, H, sa = calculate_container_properties(centers, ball_radius)
        packing_results.append({
            "balls": len(centers),
            "L": L, "W": W, "H": H,
            "sa": sa
        })
        
        # Next layer is 3.5 cm higher
        current_z += 3.5

    # Print the findings
    found_solution = False
    best_solution = None
    min_sa = initial_sa

    for res in packing_results:
        balls = res['balls']
        L, W, H, sa = res['L'], res['W'], res['H'], res['sa']
        print(f"Packing {balls} balls requires a {L:.1f}x{W:.1f}x{H:.1f} box, Surface Area: {sa:.2f} cm^2")
        
        # Check if this is a valid candidate solution
        if balls >= 27 and sa < initial_sa:
            if sa < min_sa:
                min_sa = sa
                best_solution = res
            found_solution = True

    print("-" * 30)
    
    if found_solution:
        # This part is unlikely to be reached based on the analysis
        L, W, H = best_solution['L'], best_solution['W'], best_solution['H']
        final_answer = f"{best_solution['sa']:.1f}[box {L:.1f}x{W:.1f}x{H:.1f}]"
        print("Conclusion: A more efficient container was found.")
    else:
        # The packing with 26 balls has SA=780.00, which is an improvement but doesn't hold enough balls.
        # The first packing that holds >=27 balls is the 35-ball one, but its SA (920.00) is larger.
        res_26 = packing_results[3]
        res_35 = packing_results[4]
        print(f"A box for {res_26['balls']} balls can have SA={res_26['sa']:.2f}, which is better but insufficient capacity.")
        print(f"The next step to hold >=27 balls ({res_35['balls']} balls) requires SA={res_35['sa']:.2f}, which is not an improvement.")
        print("Conclusion: No more efficient container design was found with this strategy.")
        final_answer = "0"

    print(f"\nFinal Answer: {final_answer}")

run_analysis()
<<<0>>>