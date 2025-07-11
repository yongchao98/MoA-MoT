def solve_cube_problem():
    """
    Calculates the smallest and largest possible number of green cubes.
    """
    min_total_red = float('inf')
    max_total_red = float('-inf')

    # R_e: number of red edge cubes (0 to 12)
    # R_f: number of red face-center cubes (0 to 6)
    for R_e in range(13):
        for R_f in range(7):
            # From the equation 3*R_c + 2*R_e + R_f = 18
            # We can find R_c. It must be an integer and non-negative.
            val = 18 - (2 * R_e + R_f)
            if val >= 0 and val % 3 == 0:
                R_c = val // 3
                # The number of red corners cannot exceed the total number of corners.
                if R_c <= 8:
                    
                    # A key constraint: If all 6 face-centers are red (R_f = 6),
                    # it logically forces all 12 edge cubes to be green (R_e = 0).
                    # Any combination that violates this is impossible.
                    if R_f == 6 and R_e != 0:
                        continue
                    
                    # The core cube can be red (R_i = 1) or green (R_i = 0).
                    # To find the overall min N_R, we assume the core is green (R_i = 0).
                    # To find the overall max N_R, we assume the core is red (R_i = 1).
                    
                    # Case 1: Core cube is green (R_i = 0)
                    N_R_i0 = R_c + R_e + R_f + 0
                    if N_R_i0 < min_total_red:
                        min_total_red = N_R_i0
                    if N_R_i0 > max_total_red:
                        max_total_red = N_R_i0

                    # Case 2: Core cube is red (R_i = 1)
                    N_R_i1 = R_c + R_e + R_f + 1
                    if N_R_i1 < min_total_red:
                        min_total_red = N_R_i1
                    if N_R_i1 > max_total_red:
                        max_total_red = N_R_i1

    min_green_cubes = 27 - max_total_red
    max_green_cubes = 27 - min_total_red

    print(f"The equation for red cubes is: 3 * Rc + 2 * Re + 1 * Rf = 18")
    print(f"From all valid combinations, the minimum total number of red cubes (N_R) is: {min_total_red}")
    print(f"This leads to the largest number of green cubes: 27 - {min_total_red} = {max_green_cubes}")
    print(f"The maximum total number of red cubes (N_R) is: {max_total_red}")
    print(f"This leads to the smallest number of green cubes: 27 - {max_total_red} = {min_green_cubes}")

solve_cube_problem()
<<<16 and 20>>>