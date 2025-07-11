import math

def solve_prism_permutations():
    """
    Finds the number of ordered triples (x, y, z) of positive integers
    such that xyz = 216 and x^(log2(y)) = y^(log2(z)).
    """
    V = 216
    solutions = set()

    # The condition x^(log2(y)) = y^(log2(z)) simplifies to y=1 or x=z.
    # We will find solutions for each case and use a set to store unique solutions.

    # Case 1: y = 1
    # x * 1 * z = V => x * z = V
    # x must be a divisor of V.
    count_y_eq_1 = 0
    for x in range(1, V + 1):
        if V % x == 0:
            z = V // x
            solutions.add((x, 1, z))
            count_y_eq_1 += 1
            
    # Case 2: x = z
    # x * y * x = V => x^2 * y = V
    # x^2 must be a divisor of V. We can iterate x and check.
    count_x_eq_z = 0
    # temp_solutions_x_eq_z is used to demonstrate counting for the explanation.
    temp_solutions_x_eq_z = set() 
    # Iterate x up to sqrt(V)
    for x in range(1, int(math.sqrt(V)) + 1):
        if (V % (x * x)) == 0:
            y = V // (x * x)
            solution_tuple = (x, y, x)
            solutions.add(solution_tuple)
            temp_solutions_x_eq_z.add(solution_tuple)
    count_x_eq_z = len(temp_solutions_x_eq_z)

    # The size of the set gives the total number of unique solutions.
    total_solutions = len(solutions)

    # We know the intersection is 0, but can show it programmatically
    # Intersection has y=1 and x=z => x^2 = 216, which has no integer solutions for x.
    count_intersection = count_y_eq_1 + count_x_eq_z - total_solutions
    
    print("The simplified conditions are y=1 or x=z.")
    print(f"Number of solutions where y=1: {count_y_eq_1}")
    print(f"Number of solutions where x=z: {count_x_eq_z}")
    print(f"Number of solutions satisfying both conditions (the overlap): {count_intersection}")
    print("\nUsing the Principle of Inclusion-Exclusion, the total number of permutations is:")
    print(f"{count_y_eq_1} + {count_x_eq_z} - {count_intersection} = {total_solutions}")

solve_prism_permutations()
