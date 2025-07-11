def solve_graph_problem():
    """
    This function solves the problem by using the derived equation
    2 * (b_4 - w_4) = 3 * (w_{3,R} - b_{3,R}).

    It finds the smallest positive integer for b_4 - w_4.
    """
    
    # Let D = b_4 - w_4. We are given D > 0.
    # Let K = w_{3,R} - b_{3,R}. K is an integer.
    # The equation is 2 * D = 3 * K.
    # This means 2 * D must be divisible by 3. Since 2 and 3 are coprime,
    # D must be divisible by 3.

    smallest_D = -1
    d = 1
    while True:
        # We are looking for the smallest positive integer D
        # such that 2 * D is a multiple of 3.
        if (2 * d) % 3 == 0:
            smallest_D = d
            break
        d += 1
    
    print(f"The analysis leads to the equation: 2 * (b_4 - w_4) = 3 * (w_{3,R} - b_{3,R})")
    print(f"Let D = b_4 - w_4. Since b_4 > w_4, D must be a positive integer.")
    print(f"From the equation, D must be a multiple of 3.")
    print(f"The smallest positive multiple of 3 is {smallest_D}.")
    print(f"Therefore, the smallest possible value of b_4 - w_4 is {smallest_D}.")
    print("-" * 20)
    print("To demonstrate, we can choose integer values that satisfy the condition.")
    
    # We need 2 * smallest_D = 3 * K.
    # 2 * 3 = 3 * K => K = 2.
    # So we need w_{3,R} - b_{3,R} = 2.
    # Let's pick a simple case: b_4 = 3, w_4 = 0, w_{3,R} = 2, b_{3,R} = 0.
    b4_val = smallest_D
    w4_val = 0
    w3R_val = 2
    b3R_val = 0
    
    print(f"Let's choose b_4 = {b4_val}, w_4 = {w4_val}, w_{3,R} = {w3R_val}, b_{3,R} = {b3R_val}.")
    print("The final equation with these numbers is:")
    
    left_side_val = 2 * (b4_val - w4_val)
    right_side_val = 3 * (w3R_val - b3R_val)
    
    print(f"2 * ({b4_val} - {w4_val}) = 3 * ({w3R_val} - {b3R_val})")
    print(f"2 * ({b4_val - w4_val}) = 3 * ({w3R_val - b3R_val})")
    print(f"{left_side_val} = {right_side_val}")

solve_graph_problem()