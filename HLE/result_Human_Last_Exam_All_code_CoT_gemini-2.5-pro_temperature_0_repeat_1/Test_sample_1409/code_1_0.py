import math

def solve_prism_permutations():
    """
    Calculates the number of permutations of dimensions (x, y, z) for a rectangular prism
    with volume 216, satisfying the condition x^(log_2 y) = y^(log_2 z).
    """
    V = 216
    print(f"The volume of the prism is x * y * z = {V}.")
    print("The condition x^(log_2 y) = y^(log_2 z) simplifies to (log_2 y) * (log_2 x - log_2 z) = 0.")
    print("This means a solution (x, y, z) must satisfy either y = 1 or x = z.\n")
    print("We use the Principle of Inclusion-Exclusion to find the total number of solutions.")
    print("Total = (Solutions for y=1) + (Solutions for x=z) - (Solutions for both).\n")

    # Case 1: y = 1
    # The equation is x * 1 * z = 216, so x*z = 216.
    # The number of solutions is the number of divisors of 216.
    divisors_of_V = [i for i in range(1, V + 1) if V % i == 0]
    count_y_eq_1 = len(divisors_of_V)
    print(f"1. Counting solutions where y = 1:")
    print(f"   The equation becomes x * z = {V}.")
    print(f"   The number of solutions equals the number of divisors of {V}, which is {count_y_eq_1}.")
    # print(f"   The divisors are: {divisors_of_V}")

    # Case 2: x = z
    # The equation is x * y * x = 216, so x^2 * y = 216.
    # x^2 must be a square divisor of 216.
    count_x_eq_z = 0
    solutions_x_eq_z = []
    for x in range(1, int(math.sqrt(V)) + 1):
        if V % (x * x) == 0:
            y = V // (x * x)
            solutions_x_eq_z.append(f"(x={x}, y={y}, z={x})")
            count_x_eq_z += 1
    print(f"\n2. Counting solutions where x = z:")
    print(f"   The equation becomes x^2 * y = {V}.")
    print(f"   We need to find integers x where x^2 divides {V}.")
    print(f"   The {count_x_eq_z} solutions are for x = {[1, 2, 3, 6]}.")


    # Case 3: Intersection (y = 1 and x = z)
    # The equation is x * 1 * x = 216, so x^2 = 216.
    # This requires V to be a perfect square.
    count_intersection = 0
    sqrt_V = math.isqrt(V)
    if sqrt_V * sqrt_V == V:
        count_intersection = 1
    print(f"\n3. Counting solutions in the intersection (y = 1 and x = z):")
    print(f"   The equation becomes x^2 = {V}.")
    if count_intersection == 0:
        print(f"   Since {V} is not a perfect square, there are no integer solutions for x.")
        print(f"   The number of solutions in the intersection is {count_intersection}.")
    else:
        print(f"   The number of solutions in the intersection is {count_intersection}.")


    # Final calculation using Principle of Inclusion-Exclusion
    total_permutations = count_y_eq_1 + count_x_eq_z - count_intersection
    
    print("\nFinal Calculation:")
    print(f"Total number of permutations = (Solutions for y=1) + (Solutions for x=z) - (Intersection)")
    print(f"Total = {count_y_eq_1} + {count_x_eq_z} - {count_intersection} = {total_permutations}")

solve_prism_permutations()
<<<20>>>