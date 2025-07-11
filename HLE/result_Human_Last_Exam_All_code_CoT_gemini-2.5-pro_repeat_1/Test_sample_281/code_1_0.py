def describe_drawing():
    """
    This function defines and prints the line segments of the drawing
    based on the step-by-step instructions.
    Each line is represented as a pair of (x, y) coordinates.
    """

    # Define the line segments that form the final drawing
    lines = [
        # Main square body (with right side modified by erasure)
        {'start': (0, 0), 'end': (3, 0), 'desc': "Top edge of square"},
        {'start': (0, 0), 'end': (0, -3), 'desc': "Left edge of square"},
        {'start': (0, -3), 'end': (3, -3), 'desc': "Bottom edge of square"},
        {'start': (3, 0), 'end': (3, -1), 'desc': "Upper part of right edge"},
        # Note: The segment from (3, -1) to (3, -3) is erased.

        # Pointed bottom section (like the tip of a shield)
        {'start': (0, -3), 'end': (2, -4), 'desc': "Left side of pointed bottom"},
        {'start': (2, -4), 'end': (3, -3), 'desc': "Right side of pointed bottom"},

        # Internal structure (like straps or a handle)
        # This forms a 'C' shape on the right side of the shield's interior
        {'start': (3, -1), 'end': (2, -1), 'desc': "Internal top horizontal line"},
        {'start': (2, -1), 'end': (2, -3), 'desc': "Internal vertical line"},
        {'start': (2, -3), 'end': (3, -3), 'desc': "Internal bottom horizontal line"},

        # Lines pointing to the center of the shield
        {'start': (2, -1), 'end': (1.5, -1.5), 'desc': "Line to center 1"},
        {'start': (2, -3), 'end': (1.5, -1.5), 'desc': "Line to center 2"},
    ]

    print("The final drawing consists of the following line segments:")
    for line in lines:
        start_point = line['start']
        end_point = line['end']
        description = line['desc']
        print(f"- {description}: from ({start_point[0]}, {start_point[1]}) to ({end_point[0]}, {end_point[1]})")

describe_drawing()