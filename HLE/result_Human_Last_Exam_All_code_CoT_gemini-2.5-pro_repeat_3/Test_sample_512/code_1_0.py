import math

def find_efficient_container():
    """
    Searches for a cuboid container with less surface area than the original
    12x12x12 cube that can still hold 27 spherical energy balls.
    
    The method is to pack the 27 balls in 3 layers of 9. The layers are
    then shifted horizontally relative to each other to allow for a smaller
    vertical separation, creating a denser packing.
    """
    initial_sa = 6 * 12**2
    ball_radius = 2.0
    ball_diameter = 4.0
    
    min_sa = initial_sa
    best_config = None

    # A 3x3 grid of balls with diameter 4 requires a 12x12 base.
    base_dim = 3 * ball_diameter

    # Iterate through possible horizontal layer shifts (dx, dy) on a 0.5cm grid.
    # We only need to check dx >= dy >= 0 to avoid symmetric duplicates.
    # A shift > 4cm is not useful as it allows no vertical compression.
    for dx_half in range(int(ball_diameter / 0.5) * 2): # dx from 0 to 7.5
        for dy_half in range(dx_half + 1): # dy from 0 to dx
            dx = dx_half * 0.5
            dy = dy_half * 0.5

            d_horizontal_sq = dx**2 + dy**2
            
            # If the horizontal shift is >= ball diameter, there is no overlap to exploit
            # for vertical compression. The closest balls are in different layers but
            # vertically aligned with balls in their own layer.
            if d_horizontal_sq >= ball_diameter**2:
                continue

            # Calculate the required vertical separation (dz) based on the horizontal shift.
            # The squared distance between centers of two balls in adjacent, shifted
            # layers must be >= diameter^2.
            # dx^2 + dy^2 + dz^2 >= diameter^2
            d_vertical_sq = ball_diameter**2 - d_horizontal_sq
            d_vertical = math.sqrt(d_vertical_sq)
            
            # The separation must be a multiple of 0.5 cm. We round up.
            z_sep = math.ceil(d_vertical * 2) / 2.0
            
            # If the required vertical separation is not less than the ball diameter,
            # this packing is not an improvement over the simple cubic grid.
            if z_sep >= ball_diameter:
                continue

            # We found a denser packing. Calculate the bounding box for 3 layers.
            # Layers are at z=r, z=r+z_sep, z=r+2*z_sep.
            # The first and third layers are aligned, the middle one is shifted.
            box_w = base_dim + dx
            box_h = base_dim + dy
            box_d = (2 * z_sep) + ball_diameter

            # Calculate the surface area of this new container.
            current_sa = 2 * (box_w * box_h + box_w * box_d + box_h * box_d)

            # If this is the best solution so far, save it.
            if current_sa < min_sa:
                min_sa = current_sa
                best_config = {
                    "w": box_w,
                    "h": box_h,
                    "d": box_d,
                    "sa": current_sa,
                }

    # Output the result.
    if best_config:
        sa = best_config['sa']
        w = best_config['w']
        # To maintain a consistent orientation (W>=H), we swap if needed.
        h = best_config['h'] if w >= best_config['h'] else best_config['w']
        w = best_config['w'] if w >= best_config['h'] else best_config['h']
        d = best_config['d']
        
        # Print the final equation as requested.
        print(f"SA = 2 * ({w:.1f} * {h:.1f} + {w:.1f} * {d:.1f} + {h:.1f} * {d:.1f}) = {sa:.1f}")

        # Print the final answer in the specified format.
        answer_string = f"{sa:.1f}[box {w:.1f}x{h:.1f}x{d:.1f}]"
        print(f"<<<{answer_string}>>>")
    else:
        # If no more efficient container was found.
        print("0")
        print("<<<0>>>")

# Execute the function to find the solution.
find_efficient_container()