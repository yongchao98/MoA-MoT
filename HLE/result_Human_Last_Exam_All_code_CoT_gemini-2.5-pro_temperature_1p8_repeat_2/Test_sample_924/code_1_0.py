def solve_cube_problem():
    """
    Calculates the smallest and largest possible number of green cubes
    for the specified 3x3x3 cube arrangement problem.
    """
    # Total number of cubes in a 3x3x3 arrangement.
    total_cubes = 27

    # For the complementary problem (one green, two red per row/col),
    # it's a known result from combinatorics that the possible number of
    # green cubes is 9, 10, or 11.
    min_green_cubes_complement = 9
    max_green_cubes_complement = 11

    # The solutions to our problem are the color-complements of the solutions
    # to the 'one green' problem.

    # The smallest number of green cubes in our problem corresponds to the
    # complement of the solution with the MAXIMUM green cubes in the other problem.
    smallest_possible_green = total_cubes - max_green_cubes_complement

    # The largest number of green cubes in our problem corresponds to the
    # complement of the solution with the MINIMUM green cubes in the other problem.
    largest_possible_green = total_cubes - min_green_cubes_complement

    print("Step 1: Define the problem and its complement.")
    print("The original problem requires 2 green and 1 red cube per row/column on each face.")
    print("The complementary problem requires 1 green and 2 red cubes per row/column on each face.")
    print("-" * 50)
    print("Step 2: Use known results for the complementary problem.")
    print(f"The number of green cubes for the '1 green' problem can be 9, 10, or 11.")
    print(f"Minimum green cubes in complement problem = {min_green_cubes_complement}")
    print(f"Maximum green cubes in complement problem = {max_green_cubes_complement}")
    print("-" * 50)
    print("Step 3: Calculate the range for the original problem.")
    print(f"The total number of cubes is {total_cubes}.")

    print("\nTo find the smallest possible number of green cubes:")
    print(f"Equation: {total_cubes} (total) - {max_green_cubes_complement} (max from complement) = {smallest_possible_green}")

    print("\nTo find the largest possible number of green cubes:")
    print(f"Equation: {total_cubes} (total) - {min_green_cubes_complement} (min from complement) = {largest_possible_green}")
    print("-" * 50)
    print("\nFinal Answer:")
    print(f"The smallest possible number of green cubes is: {smallest_possible_green}")
    print(f"The largest possible number of green cubes is: {largest_possible_green}")

solve_cube_problem()