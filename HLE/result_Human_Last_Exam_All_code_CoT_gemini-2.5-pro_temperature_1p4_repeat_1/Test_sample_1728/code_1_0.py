def solve_chromatic_number():
    """
    Calculates the maximal chromatic number of a graph G that is the
    sum of three cycles of length n.
    """

    # The problem asks for the *maximal* chromatic number.
    # The chromatic number of a cycle Cn is 2 if n is even, and 3 if n is odd (for n>=3).
    # To maximize the total chromatic number, we choose n to be odd for each cycle.
    max_chromatic_number_of_cycle = 3

    # The graph G is the sum of three cycles.
    num_cycles = 3

    # The chromatic number of a sum of graphs is the sum of their individual chromatic numbers.
    # So, X(G) = X(C_n) + X(C_n) + X(C_n).
    # To maximize this, we take the maximum X(C_n) for each term.
    total_max_chromatic_number = max_chromatic_number_of_cycle * num_cycles

    # We print the final equation showing each part of the sum.
    # Since G is the sum of three cycles, and the maximal chromatic number for each is 3,
    # the equation is 3 + 3 + 3.
    c1 = max_chromatic_number_of_cycle
    c2 = max_chromatic_number_of_cycle
    c3 = max_chromatic_number_of_cycle
    
    print("The maximal chromatic number is achieved when each cycle's length 'n' is odd.")
    print("The chromatic number of an odd cycle is 3.")
    print("The chromatic number of the sum of the three cycles is the sum of their individual chromatic numbers.")
    print(f"The calculation is: {c1} + {c2} + {c3} = {total_max_chromatic_number}")

solve_chromatic_number()
