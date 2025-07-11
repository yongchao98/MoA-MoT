def find_min_max_green_cubes():
    """
    Calculates the smallest and largest possible number of green cubes based on the derived constraints.
    """
    
    min_S = float('inf')
    max_S = float('-inf')

    # There are 8 corner cubes and 12 edge cubes
    # N_c is the number of green corner cubes
    # N_e is the number of green edge cubes
    for N_c in range(9):  # 0 to 8
        for N_e in range(13):  # 0 to 12
            
            # Constraint 1: Derived from the number of green face-center cubes (N_f)
            # 0 <= N_f <= 6  =>  0 <= 36 - 3*N_c - 2*N_e <= 6
            # This simplifies to: 30 <= 3*N_c + 2*N_e <= 36
            val = 3 * N_c + 2 * N_e
            if 30 <= val <= 36:
                
                S = N_c + N_e
                
                # Constraint 2: Derived from the formula for total green cubes (N_G)
                # N_G = 24 - S/2 + C_center
                # For N_G to be an integer, S must be an even number.
                if S % 2 == 0:
                    if S < min_S:
                        min_S = S
                    if S > max_S:
                        max_S = S

    # Calculate min and max N_G
    # Min N_G: Set C_center = 0 (Red) and use the largest possible S
    min_N_G = 24 - max_S / 2 + 0
    
    # Max N_G: Set C_center = 1 (Green) and use the smallest possible S
    max_N_G = 24 - min_S / 2 + 1

    print("To find the smallest number of green cubes:")
    print("We maximize the sum S = N_c + N_e, which is {}.".format(int(max_S)))
    print("We set the central cube to Red (0).")
    print("Smallest N_G = 24 - (N_c + N_e)/2 + C_center")
    print("Smallest N_G = 24 - {}/2 + 0 = {}".format(int(max_S), int(min_N_G)))
    print("-" * 20)
    print("To find the largest number of green cubes:")
    print("We minimize the sum S = N_c + N_e, which is {}.".format(int(min_S)))
    print("We set the central cube to Green (1).")
    print("Largest N_G = 24 - (N_c + N_e)/2 + C_center")
    print("Largest N_G = 24 - {}/2 + 1 = {}".format(int(min_S), int(max_N_G)))
    print("-" * 20)
    print(f"The smallest possible number of green cubes is: {int(min_N_G)}")
    print(f"The largest possible number of green cubes is: {int(max_N_G)}")


find_min_max_green_cubes()
<<<16, 19>>>