def solve_ita_matrix_code():
    """
    This script determines and prints the correct ITA Matrix extension code
    for finding direct Delta flights eligible for Platinum Companion Certificates,
    optimized for length and lexicographical order.
    """

    # Component 1: Specify the marketing airline as Delta (DL).
    # 'c dl' is shorter than 'airline dl'.
    airline_code = "c dl"

    # Component 2: Specify a direct (non-stop) flight.
    # 'n 0' is the alias and shorter than 'maxconnect 0'.
    stops_code = "n 0"

    # Component 3: Specify the eligible fare classes.
    # Platinum Companion Certificates are valid for L, U, T, X, and V fares.
    # The 'f bc=' syntax restricts the search to these booking classes per leg.
    fare_class_code = "f bc=L|U|T|X|V"

    # The problem requires minimizing string length and then maximizing the
    # case-insensitive lexicographical value. The length of the combined string
    # is the same regardless of the order of components.
    # We must arrange the components to form the lexicographically largest string.
    #
    # Comparing the components: 'n 0' > 'f bc=...' > 'c dl'.
    # Therefore, the highest lexicographical order is to start with 'n 0',
    # followed by 'f bc=...', and ending with 'c dl'.
    components = [stops_code, fare_class_code, airline_code]
    
    # The final code is the components joined by a semicolon.
    final_code = ";".join(components)

    print("The ITA Matrix outbound extension code is:")
    print(final_code)

solve_ita_matrix_code()
<<<n 0;f bc=L|U|T|X|V;c dl>>>