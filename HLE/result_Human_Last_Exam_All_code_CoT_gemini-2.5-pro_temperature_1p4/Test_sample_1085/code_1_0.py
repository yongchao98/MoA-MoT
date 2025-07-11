def solve_prime_paths():
    """
    Calculates the number of distinct Prime Paths from (1, 1) to (5, 7)
    with exactly 4 moves, interpreting a "move" as a line segment.
    """

    # Primes between the x-coordinates (1 and 5)
    primes_x = [2, 3]
    # Primes between the y-coordinates (1 and 7)
    primes_y = [2, 3, 5]

    num_choices_x = len(primes_x)
    num_choices_y = len(primes_y)

    # Calculate paths for the Horizontal-Vertical-Horizontal-Vertical (HVHV) case
    hvhv_paths = num_choices_x * num_choices_y

    # Calculate paths for the Vertical-Horizontal-Vertical-Horizontal (VHVH) case
    vhvh_paths = num_choices_y * num_choices_x

    # Total number of distinct paths is the sum of both cases
    total_paths = hvhv_paths + vhvh_paths

    # Print the explanation and the final equation
    print("This problem can be solved by interpreting a 'move' as a continuous straight-line segment.")
    print("There are two types of 4-move paths: HVHV and VHVH.")
    print(f"Number of choices for the intermediate X prime ({primes_x}): {num_choices_x}")
    print(f"Number of choices for the intermediate Y prime ({primes_y}): {num_choices_y}")
    print("\nCalculating the number of paths:")
    print(f"Number of HVHV paths = (Choices for X) * (Choices for Y) = {num_choices_x} * {num_choices_y} = {hvhv_paths}")
    print(f"Number of VHVH paths = (Choices for Y) * (Choices for X) = {num_choices_y} * {num_choices_x} = {vhvh_paths}")
    print("\nFinal equation for the total number of distinct paths:")
    # The final print statement requested by the prompt
    print(f"{hvhv_paths} + {vhvh_paths} = {total_paths}")

solve_prime_paths()
<<<12>>>