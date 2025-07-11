def solve_visibility_puzzle():
    """
    Calculates which balls Joe can see by modeling the room and obstacles
    and checking the line of sight for each ball.
    """
    print("Analyzing visibility for each ball based on room geometry.")
    print("Joe's eye can move along the line y=0.25, from x=4.5 to x=7.5, at a height z=5.\n")
    
    visible_balls = []
    
    # --- 1. Red Ball (SE Corner) ---
    print("--- 1. Red Ball (SE corner at approx. (11.75, 0.25)) ---")
    joe_x, joe_y = 7.5, 0.25
    ball_x, ball_y = 11.75, 0.25
    print(f"Considering a line of sight from Joe's rightmost position ({joe_x}, {joe_y}) to the ball at ({ball_x}, {ball_y}).")
    print("This line of sight is a straight horizontal line along y=0.25.")
    print("The wardrobe and its doors are located at y-coordinates of 4 or greater.")
    print(f"Since the sightline's y-coordinate ({joe_y}) is less than the wardrobe's minimum y-coordinate (4), the path is clear.")
    print("Result: Red ball is VISIBLE.\n")
    visible_balls.append("Red")

    # --- 2. Blue Ball (NE Corner) ---
    print("--- 2. Blue Ball (NE corner at approx. (11.75, 11.75)) ---")
    joe_x, joe_y = 4.5, 0.25  # Leftmost position to see around wardrobe
    ball_x, ball_y = 11.75, 11.75
    print(f"To see around the wardrobe, consider Joe's leftmost position ({joe_x}, {joe_y}).")
    print(f"The line of sight goes from ({joe_x}, {joe_y}) to the ball at ({ball_x}, {ball_y}).")
    print("The wardrobe's front is at x=9.5 and its north side (and north door) is at y=8.")
    
    # Equation of the line: y = m*(x - x1) + y1
    m = (ball_y - joe_y) / (ball_x - joe_x)
    y_at_wardrobe_front = m * (9.5 - joe_x) + joe_y
    
    print(f"The slope of the line is ({ball_y} - {joe_y}) / ({ball_x} - {joe_x}) = {m:.4f}.")
    print(f"Using the point-slope form: y - {joe_y} = {m:.4f} * (x - {joe_x}).")
    print(f"At the wardrobe's front (x=9.5), the sightline's y-coordinate is {m:.4f} * (9.5 - {joe_x}) + {joe_y} = {y_at_wardrobe_front:.2f}.")
    print(f"Since the sightline passes at y={y_at_wardrobe_front:.2f}, which is higher than the wardrobe's top edge at y=8, the path is clear of the wardrobe and its doors.")
    print("Result: Blue ball is VISIBLE.\n")
    visible_balls.append("Blue")

    # --- 3. Yellow Ball (NW Corner) ---
    print("--- 3. Yellow Ball (NW corner at approx. (0.25, 11.75)) ---")
    joe_x, joe_y = 7.5, 0.25  # Rightmost position to see around bookshelf
    ball_x, ball_y = 0.25, 11.75
    print(f"To see around the bookshelf, consider Joe's rightmost position ({joe_x}, {joe_y}).")
    print(f"The line of sight goes from ({joe_x}, {joe_y}) to the ball at ({ball_x}, {ball_y}).")
    print("The bookshelf occupies the area from x=0 to x=1 and y=0 to y=4.")
    
    m = (ball_y - joe_y) / (ball_x - joe_x)
    y_at_bookshelf_edge = m * (1.0 - joe_x) + joe_y

    print(f"The slope of the line is ({ball_y} - {joe_y}) / ({ball_x} - {joe_x}) = {m:.4f}.")
    print(f"Using the point-slope form: y - {joe_y} = {m:.4f} * (x - {joe_x}).")
    print(f"At the bookshelf's front edge (x=1.0), the sightline's y-coordinate is {m:.4f} * (1.0 - {joe_x}) + {joe_y} = {y_at_bookshelf_edge:.2f}.")
    print(f"Since this y-value ({y_at_bookshelf_edge:.2f}) is far north of the bookshelf's maximum extent (y=4), the path is clear.")
    print("Result: Yellow ball is VISIBLE.\n")
    visible_balls.append("Yellow")

    # --- 4. Green Ball (on SW shelf) ---
    print("--- 4. Green Ball (SW corner on shelf, approx. (0.25, 0.25)) ---")
    ball_x = 0.25
    door_x = 4.5
    print(f"The ball is at an x-coordinate of {ball_x}.")
    print(f"The room door is hinged at x={door_x}, forming a solid obstacle taller than Joe.")
    print(f"Joe stands in the doorway at an x-coordinate greater than {door_x} (e.g., 4.51 to 7.5).")
    print(f"Any line of sight from Joe (x > {door_x}) to the ball (x < {door_x}) must pass through the plane of the door at x={door_x}, which blocks the view.")
    print("Result: Green ball is NOT VISIBLE.\n")

    # --- 5. Purple Ball (in wardrobe) ---
    print("--- 5. Purple Ball (in SE corner of wardrobe, approx. (11.75, 4.25)) ---")
    joe_x, joe_y = 4.5, 0.25  # Most extreme viewpoint to peer inside
    ball_x, ball_y = 11.75, 4.25
    print(f"Consider Joe's leftmost position ({joe_x}, {joe_y}) to get the best angle into the wardrobe.")
    print("The wardrobe opening is at x=9.5, from y=4 to y=8.")
    
    m = (ball_y - joe_y) / (ball_x - joe_x)
    y_at_wardrobe_front = m * (9.5 - joe_x) + joe_y

    print(f"The slope of the line is ({ball_y} - {joe_y}) / ({ball_x} - {joe_x}) = {m:.4f}.")
    print(f"Using the point-slope form: y - {joe_y} = {m:.4f} * (x - {joe_x}).")
    print(f"At the wardrobe's front (x=9.5), the sightline's y-coordinate is {m:.4f} * (9.5 - {joe_x}) + {joe_y} = {y_at_wardrobe_front:.2f}.")
    print(f"This ({y_at_wardrobe_front:.2f}) is below the wardrobe opening which starts at y=4. The line of sight is blocked by the solid front of the wardrobe.")
    print("Result: Purple ball is NOT VISIBLE.\n")
    
    final_answer = ", ".join(visible_balls)
    print(f"Final conclusion: Joe can see the {final_answer} balls.")
    
    print(f'<<<{final_answer}>>>')

solve_visibility_puzzle()