import matplotlib.pyplot as plt

def draw_figure():
    """
    This function draws the figure based on the textual description.
    The interpretation of the steps is explained in the comments.
    """
    # Create a figure and axes for plotting
    fig, ax = plt.subplots()

    # --- Step 1-4: The main square body ---
    # The square has corners at (0,0), (3,0), (3,3), (0,3).
    # b1 is (0,3), b2 is (3,3).
    # We will draw the parts that are not erased later.
    top_wall = ([0, 3], [0, 0])
    left_wall = ([0, 0], [0, 3])
    floor = ([0, 3], [3, 3])
    # The right wall is drawn in two parts because a section is erased.
    # The erased section is from y=1 to y=3.
    right_wall_top = ([3, 3], [0, 1])
    
    ax.plot(top_wall[0], top_wall[1], 'b-')
    ax.plot(left_wall[0], left_wall[1], 'b-')
    ax.plot(floor[0], floor[1], 'b-')
    ax.plot(right_wall_top[0], right_wall_top[1], 'b-')
    
    # --- Step 5-7: The roof ---
    # b1=(0,3), b2=(3,3). The peak of the roof, p, is at (2,4).
    b1 = (0, 3)
    b2 = (3, 3)
    p = (2, 4)
    left_roof = ([b1[0], p[0]], [b1[1], p[1]])
    right_roof = ([p[0], b2[0]], [p[1], b2[1]])
    
    ax.plot(left_roof[0], left_roof[1], 'b-')
    ax.plot(right_roof[0], right_roof[1], 'b-')

    # --- Step 8-15: The side structure and internal lines ---
    # These steps have contradictions, so we follow a plausible interpretation.
    # c is the center of the square at (1.5, 1.5)
    c = (1.5, 1.5)
    # a1 is at (2,1), a2 is at (2,3)
    a1 = (2, 1)
    a2 = (2, 3)
    # s is at (3,1)
    s = (3, 1)
    
    # Draw segments connecting a1 to c and a2 to c (Step 13)
    line_a1_c = ([a1[0], c[0]], [a1[1], c[1]])
    line_a2_c = ([a2[0], c[0]], [a2[1], c[1]])
    ax.plot(line_a1_c[0], line_a1_c[1], 'b-')
    ax.plot(line_a2_c[0], line_a2_c[1], 'b-')

    # The rectangle from step 14 has corners a1, a2, b2, s.
    # The segment s-b2 is erased (Step 15), so we draw the other 3 sides.
    # This forms a U-shape attached to the house.
    rect_top = ([a1[0], s[0]], [a1[1], s[1]])
    rect_left = ([a1[0], a2[0]], [a1[1], a2[1]])
    rect_bottom = ([a2[0], b2[0]], [a2[1], b2[1]])
    
    ax.plot(rect_top[0], rect_top[1], 'b-')
    ax.plot(rect_left[0], rect_left[1], 'b-')
    ax.plot(rect_bottom[0], rect_bottom[1], 'b-')

    # --- Finalizing the plot ---
    # Set plot limits and aspect ratio
    ax.set_xlim(-1, 5)
    ax.set_ylim(5, -1)  # Invert y-axis to match "lined paper" description
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks(range(-1, 6))
    ax.set_yticks(range(-1, 6))
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.set_title("Resulting Drawing")
    
    # Print the coordinates of the key points
    print("Key points coordinates (x, y):")
    print(f"b1 (bottom-left of square): {b1}")
    print(f"b2 (bottom-right of square): {b2}")
    print(f"p (roof peak): {p}")
    print(f"c (center of square): {c}")
    print(f"a1 (top-left of side structure): {a1}")
    print(f"a2 (bottom-left of side structure): {a2}")
    
    plt.show()

# Execute the function to generate the plot
draw_figure()