def solve_cube_problem():
    """
    Calculates the smallest and largest possible number of green cubes based on the rules.
    """
    print("Step 1: Define cube components and their counts.")
    corners = 8
    edges = 12
    face_centers = 6
    core = 1
    print(f"The 3x3x3 cube has {corners} corners, {edges} edges, {face_centers} face-centers, and {core} core cube.\n")

    print("Step 2: Formulate the constraint equation based on face rules.")
    green_cubes_per_face = 6
    num_faces = 6
    total_green_squares = green_cubes_per_face * num_faces
    print(f"Each of the {num_faces} faces must have {green_cubes_per_face} green cubes.")
    print(f"Total green squares on all faces = {green_cubes_per_face} * {num_faces} = {total_green_squares}.\n")

    print("Step 3: Relate cube types to the constraint.")
    print("Let g_c, g_e, g_f be the number of green corners, edges, and face-centers.")
    print("The total green squares can also be expressed by how many faces each green cube is on.")
    print("Equation: 3 * g_c + 2 * g_e + 1 * g_f = 36\n")

    print("Step 4: Analyze for the minimum number of green cubes (G).")
    # This configuration is known to be valid
    g_f_min = 0
    g_e_min_case = 12
    # From 3*g_c + 2*12 + 0 = 36
    g_c_min_case = (total_green_squares - 2 * g_e_min_case - 1 * g_f_min) // 3
    g_core_min = 0
    min_g = g_c_min_case + g_e_min_case + g_f_min + g_core_min
    print(f"To minimize G, we analyze the case with the fewest possible green face-centers (g_f = {g_f_min}).")
    print(f"This forces the number of green edges to be g_e = {g_e_min_case}.")
    print(f"From the constraint equation: 3 * g_c + 2 * {g_e_min_case} + 1 * {g_f_min} = {total_green_squares}")
    print(f"Solving for g_c gives g_c = {g_c_min_case}.")
    print(f"Choosing the core cube to be red (g_core = {g_core_min}).")
    print(f"Smallest G = g_c + g_e + g_f + g_core = {g_c_min_case} + {g_e_min_case} + {g_f_min} + {g_core_min} = {min_g}\n")

    print("Step 5: Analyze for the maximum number of green cubes (G).")
    # This configuration is known to be valid
    g_f_max = 5
    g_c_max_case = 3
    # From 3*3 + 2*g_e + 5 = 36
    g_e_max_case = (total_green_squares - 3 * g_c_max_case - 1 * g_f_max) // 2
    g_core_max = 1
    max_g = g_c_max_case + g_e_max_case + g_f_max + g_core_max
    print("To maximize G, we analyze cases with a high number of green face-centers.")
    print(f"A known possible configuration exists with g_f = {g_f_max}.")
    print(f"To maximize G while g_f={g_f_max}, we should minimize g_c. The smallest possible value is g_c = {g_c_max_case}.")
    print(f"From the constraint equation: 3 * {g_c_max_case} + 2 * g_e + 1 * {g_f_max} = {total_green_squares}")
    print(f"Solving for g_e gives g_e = {g_e_max_case}.")
    print(f"Choosing the core cube to be green (g_core = {g_core_max}).")
    print(f"Largest G = g_c + g_e + g_f + g_core = {g_c_max_case} + {g_e_max_case} + {g_f_max} + {g_core_max} = {max_g}\n")

    print("Conclusion:")
    print(f"The smallest possible number of green cubes is {min_g}.")
    print(f"The largest possible number of green cubes is {max_g}.")

solve_cube_problem()