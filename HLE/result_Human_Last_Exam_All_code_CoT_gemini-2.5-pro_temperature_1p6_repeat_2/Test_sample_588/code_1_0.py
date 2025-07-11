def find_solution_tuple():
    """
    This function programmatically finds the lexicographically least tuple
    (a1, b1, a2, b2, a3, b3) satisfying the derived topological conditions.
    The minimal number of manifolds, l, is 3. The governing equation is:
    (1-a1)(1-b1) + (1-a2)(1-b2) + (1-a3)(1-b3) = 1
    where ai, bi are non-negative integers not equal to 1.
    """
    # Define a reasonable search space for a_i and b_i.
    # The smallest values are {0, 2, 3, ...}.
    search_limit = 5
    possible_values = [0] + list(range(2, search_limit))

    # We search for the tuple by iterating through possibilities in lexicographical order.
    # To simplify, we can assume a_i <= b_i for each pair.
    # We also enforce that the pairs are sorted: (a1,b1) <= (a2,b2) <= (a3,b3).
    print("Searching for the lexicographically least tuple (a1,b1,a2,b2,a3,b3)...")

    for a1 in possible_values:
        for b1 in [v for v in possible_values if v >= a1]:
            for a2 in possible_values:
                for b2 in [v for v in possible_values if v >= a2]:
                    # Enforce lexicographical ordering of pairs
                    if (a2, b2) < (a1, b1):
                        continue

                    for a3 in possible_values:
                        for b3 in [v for v in possible_values if v >= a3]:
                            # Enforce lexicographical ordering of pairs
                            if (a3, b3) < (a2, b2):
                                continue

                            # Calculate the terms of the equation
                            x1 = (1 - a1) * (1 - b1)
                            x2 = (1 - a2) * (1 - b2)
                            x3 = (1 - a3) * (1 - b3)

                            # Check if the equation is satisfied
                            if x1 + x2 + x3 == 1:
                                result_tuple = (a1, b1, a2, b2, a3, b3)
                                
                                print("\nFound solution:")
                                print("Equation: x1 + x2 + x3 = 1")
                                print(f"Pair 1: ({a1}, {b1}) -> x1 = (1-{a1}) * (1-{b1}) = {x1}")
                                print(f"Pair 2: ({a2}, {b2}) -> x2 = (1-{a2}) * (1-{b2}) = {x2}")
                                print(f"Pair 3: ({a3}, {b3}) -> x3 = (1-{a3}) * (1-{b3}) = {x3}")
                                print(f"Check sum: {x1} + {x2} + {x3} = {x1 + x2 + x3}")
                                
                                # Format the final answer as a flat tuple with no spaces
                                answer_string = str(result_tuple).replace(" ", "")
                                print(f"\nThe lexicographically least tuple is {answer_string}")
                                return

# Run the search function
find_solution_tuple()
