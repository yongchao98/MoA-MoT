def solve_cube_problem():
    """
    This function calculates the smallest and largest possible number of green cubes
    based on the derived mathematical constraints.
    """
    total_cubes = 27

    print("The problem is solved by setting up and solving a system of linear Diophantine equations based on the cube's geometry and rules.")
    print("The key equations derived are:")
    print("1) 3*n_R_c + n_R_e = 12")
    print("2) n_R_e + n_R_f = 6")
    print("3) N_R = n_R_c + n_R_e + n_R_f + n_R_core")
    print("where n_R_c, n_R_e, n_R_f are the numbers of red corner, edge, and face-center cubes, and N_R is the total number of red cubes.\n")

    # Find the smallest number of green cubes (n_G_min), which corresponds to the largest number of red cubes (N_R_max).
    # To maximize N_R = n_R_c + 6 + n_R_core, we need to maximize n_R_c and set n_R_core=1.
    # The valid range for n_R_c is [2, 4]. We pick the maximum.
    n_rc_max = 4
    n_re_for_max_nr = 12 - 3 * n_rc_max
    n_rf_for_max_nr = 6 - n_re_for_max_nr
    n_rcore_for_max_nr = 1 # Make core red to maximize red cubes

    N_R_max = n_rc_max + n_re_for_max_nr + n_rf_for_max_nr + n_rcore_for_max_nr
    n_G_min = total_cubes - N_R_max
    
    print("--- Smallest Number of Green Cubes ---")
    print(f"To minimize green cubes, we maximize red cubes (N_R). This happens when n_R_c is maximum.")
    print(f"Set n_R_c = {n_rc_max}.")
    print(f"This implies n_R_e = 12 - 3 * {n_rc_max} = {n_re_for_max_nr}.")
    print(f"And n_R_f = 6 - n_R_e = 6 - {n_re_for_max_nr} = {n_rf_for_max_nr}.")
    print(f"We also make the core cube red, so n_R_core = {n_rcore_for_max_nr}.")
    print(f"Maximum N_R = n_R_c + n_R_e + n_R_f + n_R_core = {n_rc_max} + {n_re_for_max_nr} + {n_rf_for_max_nr} + {n_rcore_for_max_nr} = {N_R_max}")
    print(f"Smallest number of green cubes = {total_cubes} - {N_R_max} = {n_G_min}\n")

    # Find the largest number of green cubes (n_G_max), which corresponds to the smallest number of red cubes (N_R_min).
    # To minimize N_R = n_R_c + 6 + n_R_core, we need to minimize n_R_c and set n_R_core=0.
    # The valid range for n_R_c is [2, 4]. We pick the minimum.
    n_rc_min = 2
    n_re_for_min_nr = 12 - 3 * n_rc_min
    n_rf_for_min_nr = 6 - n_re_for_min_nr
    n_rcore_for_min_nr = 0 # Make core green to minimize red cubes

    N_R_min = n_rc_min + n_re_for_min_nr + n_rf_for_min_nr + n_rcore_for_min_nr
    n_G_max = total_cubes - N_R_min
    
    print("--- Largest Number of Green Cubes ---")
    print(f"To maximize green cubes, we minimize red cubes (N_R). This happens when n_R_c is minimum.")
    print(f"Set n_R_c = {n_rc_min}.")
    print(f"This implies n_R_e = 12 - 3 * {n_rc_min} = {n_re_for_min_nr}.")
    print(f"And n_R_f = 6 - n_R_e = 6 - {n_re_for_min_nr} = {n_rf_for_min_nr}.")
    print(f"We also make the core cube green, so n_R_core = {n_rcore_for_min_nr}.")
    print(f"Minimum N_R = n_R_c + n_R_e + n_R_f + n_R_core = {n_rc_min} + {n_re_for_min_nr} + {n_rf_for_min_nr} + {n_rcore_for_min_nr} = {N_R_min}")
    print(f"Largest number of green cubes = {total_cubes} - {N_R_min} = {n_G_max}\n")

solve_cube_problem()