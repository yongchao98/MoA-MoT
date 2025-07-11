import math

def solve_prism_permutations():
    """
    This function calculates the number of permutations of dimensions (x, y, z)
    of a rectangular prism with volume 216, satisfying the condition
    x^(log_2 y) = y^(log_2 z).
    
    The condition simplifies to either y=1 or x=z.
    """
    volume = 216
    
    # We use a set to store unique solutions (ordered triplets).
    solutions = set()

    # Find all divisors of the volume
    divisors = []
    for i in range(1, volume + 1):
        if volume % i == 0:
            divisors.append(i)

    # --- Case 1: y = 1 ---
    # The volume equation becomes x * 1 * z = 216.
    # The number of solutions is the number of divisors of 216.
    solutions_y_1 = []
    for x in divisors:
        z = volume // x
        solutions.add((x, 1, z))
        solutions_y_1.append((x, 1, z))
        
    print("Solutions from the condition y = 1:")
    # Printing sorted list for better readability
    for s in sorted(solutions_y_1):
        print(s)
    num_solutions_y_1 = len(solutions_y_1)
    print(f"\nNumber of solutions where y = 1: {num_solutions_y_1}")
    print("-" * 30)

    # --- Case 2: x = z ---
    # The volume equation is x^2 * y = 216.
    # x^2 must be a perfect square that divides 216.
    solutions_x_eq_z = []
    for d in divisors:
        # Check if the divisor is a perfect square
        sqrt_d = int(math.sqrt(d))
        if sqrt_d * sqrt_d == d:
            x = sqrt_d
            y = volume // d
            solutions.add((x, y, x))
            solutions_x_eq_z.append((x, y, x))

    print("Solutions from the condition x = z:")
    # Printing sorted list for better readability
    for s in sorted(solutions_x_eq_z):
        print(s)
    num_solutions_x_eq_z = len(solutions_x_eq_z)
    print(f"\nNumber of solutions where x = z: {num_solutions_x_eq_z}")
    print("-" * 30)
    
    # The two sets of solutions are disjoint, so the total is their sum.
    total_solutions = num_solutions_y_1 + num_solutions_x_eq_z
    
    print("Final Calculation:")
    print(f"The total number of permutations is the sum of counts from both cases.")
    print(f"Total = (solutions where y = 1) + (solutions where x = z)")
    print(f"Total = {num_solutions_y_1} + {num_solutions_x_eq_z} = {total_solutions}")

solve_prism_permutations()
<<<20>>>