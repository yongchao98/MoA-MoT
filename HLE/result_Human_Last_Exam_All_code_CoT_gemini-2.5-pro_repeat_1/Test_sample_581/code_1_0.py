def solve_cap_set_lower_bound():
    """
    This script calculates a lower bound for the size of a cap set in dimension 8
    and then presents the best-known lower bound from mathematical research.

    A cap set is a collection of points in the vector space (Z/3)^n such that no
    three distinct points are collinear (i.e., for any distinct x, y, z in the set, x + y + z != 0).
    The size of the largest possible cap set in dimension n is denoted by r_3(n).

    Lower bounds can be found using the product construction rule:
    r_3(a + b) >= r_3(a) * r_3(b).
    """

    # These are the exact sizes of cap sets for dimensions n=1 to 6.
    r3_known_values = {
        1: 2,
        2: 4,
        3: 9,
        4: 20,
        5: 45,
        6: 112,
    }

    print("--- Calculating Lower Bounds for r_3(8) using Product Construction ---")
    print("The rule is: r_3(a + b) >= r_3(a) * r_3(b)\n")

    n = 8
    best_product_bound = 0
    best_equation = ""

    # We check all partitions of 8 into two numbers 'a' and 'b' for which we know r_3(a) and r_3(b).
    for a in range(1, n // 2 + 1):
        b = n - a
        if a in r3_known_values and b in r3_known_values:
            bound = r3_known_values[a] * r3_known_values[b]
            
            # The prompt asks to output each number in the final equation.
            # We will show the equation for each product calculation.
            equation_str = f"For n = {a} + {b}: r_3({a}) * r_3({b}) = {r3_known_values[a]} * {r3_known_values[b]} = {bound}"
            print(equation_str)

            if bound > best_product_bound:
                best_product_bound = bound
                best_equation = equation_str

    print(f"\nThe best lower bound found using the simple product construction is {best_product_bound}.")
    print(f"This comes from the calculation: {best_equation.split(': ')[1]}")

    print("\n--- Best Known Lower Bound ---")
    print("While the product construction gives a bound of 448, more advanced constructions")
    print("have established a better lower bound.")
    
    # This is the current record from mathematical literature (Y. Edel, 2017).
    best_known_lower_bound = 496
    
    print("\nThe best known lower bound for the size of a cap set in dimension 8 is:")
    print(best_known_lower_bound)

solve_cap_set_lower_bound()