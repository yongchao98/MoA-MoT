import math

def solve_prism_permutations():
    """
    Calculates the number of permutations of dimensions (x, y, z) for a rectangular prism
    with volume 216, satisfying the condition x^(log2 y) = y^(log2 z).
    """
    V = 216

    print(f"The volume of the prism is xyz = {V}.")
    print("The condition x^(log_2 y) = y^(log_2 z) simplifies to (y=1) or (x=z).")
    print("\nWe use the Principle of Inclusion-Exclusion to find the total number of permutations (x, y, z).")
    print("Total = (Number of solutions where y=1) + (Number of solutions where x=z) - (Number of solutions where y=1 and x=z).")

    # --- Case 1: Count solutions where y = 1 ---
    # This means x * z = 216. The number of solutions (x, z)
    # is the number of divisors of 216.
    def get_divisors_count(n):
        count = 0
        for i in range(1, int(math.sqrt(n)) + 1):
            if n % i == 0:
                if n // i == i:
                    count += 1
                else:
                    count += 2
        return count

    count_y_is_1 = get_divisors_count(V)

    print("\n--- Case 1: y = 1 ---")
    print(f"The volume equation becomes x * 1 * z = {V}, so x*z = {V}.")
    print(f"The number of integer solutions (x, z) is the number of divisors of {V}.")
    print(f"Number of solutions where y=1 is {count_y_is_1}.")

    # --- Case 2: Count solutions where x = z ---
    # This means x^2 * y = 216.
    # x^2 must be a square factor of 216.
    count_x_is_z = 0
    solutions_x_is_z = []
    # We only need to check for x up to floor(sqrt(V))
    for x in range(1, int(V**0.5) + 1):
        if (x * x) != 0 and V % (x * x) == 0:
            count_x_is_z += 1
            y = V // (x * x)
            solutions_x_is_z.append((x, y, x))
    
    print("\n--- Case 2: x = z ---")
    print(f"The volume equation becomes x * y * x = {V}, so x^2 * y = {V}.")
    print(f"This requires x^2 to be a perfect square divisor of {V}.")
    print(f"The {count_x_is_z} solutions (x, y, z) are:")
    for sol in solutions_x_is_z:
        print(f"  ({sol[0]}, {sol[1]}, {sol[2]})")
    print(f"Number of solutions where x=z is {count_x_is_z}.")


    # --- Case 3: Count solutions where y = 1 AND x = z ---
    # This means x^2 * 1 = 216.
    count_intersection = 0
    sqrt_V = math.isqrt(V)
    if sqrt_V * sqrt_V == V:
        count_intersection = 1
    
    print("\n--- Case 3: Intersection (y = 1 and x = z) ---")
    print(f"This requires x^2 * 1 = {V}, so x = sqrt({V}).")
    if count_intersection == 0:
        print(f"Since sqrt({V}) is not an integer, there are no integer solutions for x.")
        print("Number of solutions in the intersection is 0.")
    else:
        print(f"Since sqrt({V}) = {sqrt_V} is an integer, there is one solution.")
        print("Number of solutions in the intersection is 1.")

    # --- Final Calculation ---
    total_permutations = count_y_is_1 + count_x_is_z - count_intersection
    print("\n--- Final Calculation ---")
    print("Total number of permutations = (Solutions for y=1) + (Solutions for x=z) - (Solutions for both)")
    print(f"Total = {count_y_is_1} + {count_x_is_z} - {count_intersection}")
    print(f"The total number of permutations is {total_permutations}.")

solve_prism_permutations()