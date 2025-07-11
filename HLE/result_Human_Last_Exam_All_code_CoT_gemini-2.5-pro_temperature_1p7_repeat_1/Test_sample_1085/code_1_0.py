import math

def solve_prime_paths():
    """
    Calculates the number of distinct 4-move Prime Paths from (1,1) to (5,7).
    """

    # Step 1: Define the sets of available intermediate prime coordinates.
    # The destination is (5,7). Paths start at (1,1).
    # Primes up to 7 are {2, 3, 5, 7}.
    # Intermediate x-coordinates cannot be 1 or 5.
    intermediate_x_primes = {2, 3, 7}
    num_x_choices = len(intermediate_x_primes)

    # Intermediate y-coordinates cannot be 1 or 7.
    intermediate_y_primes = {2, 3, 5}
    num_y_choices = len(intermediate_y_primes)

    print("Step-by-step calculation of distinct Prime Paths from (1,1) to (5,7):\n")

    # Step 2: Case 1: 2 Horizontal (H) and 2 Vertical (V) moves.
    # This path structure requires choosing one intermediate x and one intermediate y coordinate.
    # The number of ways to arrange 2 H and 2 V moves is C(4,2).
    c_4_2 = math.comb(4, 2)
    # The number of ways to choose the intermediate x and y primes are independent.
    paths_2h_2v = c_4_2 * num_x_choices * num_y_choices
    print(f"Case 1 (2 Horizontal, 2 Vertical moves):")
    print(f"  Number of move patterns (e.g., HVHV): C(4, 2) = {c_4_2}")
    print(f"  Number of choices for intermediate x-prime ({intermediate_x_primes}): {num_x_choices}")
    print(f"  Number of choices for intermediate y-prime ({intermediate_y_primes}): {num_y_choices}")
    print(f"  Subtotal = {c_4_2} * {num_x_choices} * {num_y_choices} = {paths_2h_2v}\n")


    # Step 3: Case 2: 1 Horizontal (H) and 3 Vertical (V) moves.
    # This requires choosing two distinct intermediate y-coordinates in a specific order.
    # The number of ways to arrange 1 H and 3 V moves is C(4,1).
    c_4_1 = math.comb(4, 1)
    # The number of ways to pick an ordered pair of y-primes from 3 choices is P(3,2).
    p_3_2_y = math.perm(num_y_choices, 2)
    paths_1h_3v = c_4_1 * p_3_2_y
    print(f"Case 2 (1 Horizontal, 3 Vertical moves):")
    print(f"  Number of move patterns (e.g., HVVV): C(4, 1) = {c_4_1}")
    print(f"  Number of ordered choices for 2 intermediate y-primes from {num_y_choices}: P(3, 2) = {p_3_2_y}")
    print(f"  Subtotal = {c_4_1} * {p_3_2_y} = {paths_1h_3v}\n")


    # Step 4: Case 3: 3 Horizontal (H) and 1 Vertical (V) move.
    # This requires choosing two distinct intermediate x-coordinates in a specific order.
    # The number of ways to arrange 3 H and 1 V moves is C(4,3).
    c_4_3 = math.comb(4, 3)
    # The number of ways to pick an ordered pair of x-primes from 3 choices is P(3,2).
    p_3_2_x = math.perm(num_x_choices, 2)
    paths_3h_1v = c_4_3 * p_3_2_x
    print(f"Case 3 (3 Horizontal, 1 Vertical move):")
    print(f"  Number of move patterns (e.g., VHHH): C(4, 3) = {c_4_3}")
    print(f"  Number of ordered choices for 2 intermediate x-primes from {num_x_choices}: P(3, 2) = {p_3_2_x}")
    print(f"  Subtotal = {c_4_3} * {p_3_2_x} = {paths_3h_1v}\n")

    # Step 5: Final calculation
    total_paths = paths_2h_2v + paths_1h_3v + paths_3h_1v
    print("Final Calculation:")
    print(f"  Total Paths = (Subtotal from Case 1) + (Subtotal from Case 2) + (Subtotal from Case 3)")
    print(f"  Total Paths = {paths_2h_2v} + {paths_1h_3v} + {paths_3h_1v} = {total_paths}")


solve_prime_paths()