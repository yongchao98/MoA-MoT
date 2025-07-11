import math

def solve_prism_permutations():
    """
    Finds the number of permutations of dimensions (x, y, z) for a rectangular prism
    with volume 216, satisfying x^(log_2 y) = y^(log_2 z).
    """
    volume = 216

    print("The problem is to find the number of ordered triples (x, y, z) of positive integers such that:")
    print(f"1. x * y * z = {volume}")
    print("2. x^(log_2 y) = y^(log_2 z)")
    print("\nThe second condition simplifies to (y = 1) or (x = z).")
    print("We use the Principle of Inclusion-Exclusion to find the total number of solutions.\n")

    # --- Case 1: Count solutions where y = 1 ---
    # The volume equation is x * 1 * z = 216.
    # The number of ordered integer pairs (x, z) is the number of divisors of 216.
    count_y_is_1 = 0
    for i in range(1, volume + 1):
        if volume % i == 0:
            count_y_is_1 += 1
    
    print("Step 1: Finding solutions where y = 1.")
    print(f"The equation becomes x * z = {volume}.")
    print(f"The number of ordered pairs (x, z) is the number of divisors of {volume}.")
    print(f"Number of solutions with y = 1 is {count_y_is_1}.")
    print("-" * 30)

    # --- Case 2: Count solutions where x = z ---
    # The volume equation is x * y * x = 216, or x^2 * y = 216.
    # We need to find integer pairs (x, y) that satisfy this.
    # This means x^2 must be a divisor of 216.
    count_x_equals_z = 0
    solutions_x_equals_z = []
    # We only need to check x up to sqrt(216)
    limit = int(math.sqrt(volume))
    for x in range(1, limit + 2): # Loop a bit beyond the limit for safety
        if (x * x) > volume:
            break
        if volume % (x * x) == 0:
            y = volume // (x * x)
            solutions_x_equals_z.append((x, y, x))
            count_x_equals_z += 1
    
    print("Step 2: Finding solutions where x = z.")
    print(f"The equation becomes x^2 * y = {volume}.")
    print(f"We need to find integers x where x^2 divides {volume}.")
    print(f"The solutions (x, y, z) are: {solutions_x_equals_z}")
    print(f"Number of solutions with x = z is {count_x_equals_z}.")
    print("-" * 30)

    # --- Case 3: Count solutions where y = 1 AND x = z (Intersection) ---
    # The volume equation is x * 1 * x = 216, or x^2 = 216.
    # We check if 216 is a perfect square.
    count_intersection = 0
    sqrt_volume = math.sqrt(volume)
    if sqrt_volume == int(sqrt_volume):
        count_intersection = 1
    
    print("Step 3: Finding solutions where y = 1 AND x = z.")
    print(f"The equation becomes x^2 = {volume}.")
    print(f"Since sqrt({volume}) = {sqrt_volume:.4f}, which is not an integer, there are no integer solutions.")
    print(f"Number of solutions in the intersection is {count_intersection}.")
    print("-" * 30)

    # --- Final Calculation ---
    total_permutations = count_y_is_1 + count_x_equals_z - count_intersection
    
    print("Final Calculation using the Principle of Inclusion-Exclusion:")
    print("Total Permutations = (Solutions with y=1) + (Solutions with x=z) - (Intersection)")
    print(f"Total Permutations = {count_y_is_1} + {count_x_equals_z} - {count_intersection}")
    print(f"\nThe total number of permutations is: {total_permutations}")

solve_prism_permutations()