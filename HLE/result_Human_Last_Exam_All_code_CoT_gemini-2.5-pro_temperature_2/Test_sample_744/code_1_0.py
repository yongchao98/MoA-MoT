def solve_puzzle():
    """
    Calculates the maximum value that can be produced from a square hollow tube.

    The strategy is to make as many 'whole balls' as possible, as they provide
    a higher value for the same amount of material space compared to 'welded balls'.
    """

    # --- Given Parameters ---
    outer_side_cm = 20
    length_cm = 100
    thickness_cm = 4
    ball_radius_cm = 2
    value_whole_ball = 3

    # --- Calculations ---

    # A whole ball has a 4cm diameter and is cut from a 4x4x4 cm cube.
    ball_diameter_cm = ball_radius_cm * 2

    # Calculate the number of balls that can be packed in the tube's cross-section.
    # First, find the inner dimension of the hollow tube.
    inner_side_cm = outer_side_cm - 2 * thickness_cm
    # Then, calculate the area of the material in the cross-section.
    material_cross_section_area = outer_side_cm**2 - inner_side_cm**2
    # Then, find the footprint area of one ball's bounding cube.
    ball_footprint_area = ball_diameter_cm**2
    # The number of balls that fit side-by-side in one slice is the ratio of these areas.
    balls_per_slice = material_cross_section_area // ball_footprint_area

    # Calculate how many 4cm-thick slices of balls can be made along the tube's length.
    num_slices = length_cm // ball_diameter_cm

    # Calculate the total number of whole balls that can be manufactured.
    total_whole_balls = balls_per_slice * num_slices

    # Calculate the final maximum value.
    max_value = total_whole_balls * value_whole_ball

    # --- Output the result ---
    print("The optimal strategy is to produce only whole balls.")
    print(f"First, we find the total number of whole balls that can be made:")
    print(f"Total Balls = (Balls per Cross-Section Slice) * (Number of Slices)")
    print(f"Total Balls = ({int(balls_per_slice)}) * ({int(num_slices)}) = {int(total_whole_balls)}")
    print("\nThen, we calculate the maximum value:")
    print("Final Equation:")
    print(f"{int(total_whole_balls)} whole balls * {value_whole_ball} value/ball = {max_value}")

solve_puzzle()