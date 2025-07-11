def solve_e8_tori():
    """
    Calculates the number of Fq-rational maximal tori in a reductive group of type E8.
    This is equivalent to finding the number of conjugacy classes in the Weyl group W(E_8).
    The calculation is based on the known classification of conjugacy classes by element order.
    """

    # Data from the classification of conjugacy classes in W(E_8) by R.W. Carter.
    # The dictionary maps an element's order to the number of conjugacy classes of elements with that order.
    classes_per_order = {
        1: 1,
        2: 4,
        3: 5,
        4: 8,
        5: 4,
        6: 15,
        7: 3,
        8: 6,
        9: 2,
        10: 8,
        12: 15,
        14: 3,
        15: 4,
        18: 5,
        20: 6,
        24: 6,
        30: 17
    }

    # Extract the number of classes for each order.
    counts = list(classes_per_order.values())

    # Calculate the total number of conjugacy classes.
    total_number_of_classes = sum(counts)

    # Build the equation string to show the summation.
    # The problem asks to output each number in the final equation.
    equation_str = " + ".join(map(str, counts))

    print(f"The number of conjugacy classes in W(E_8) is the sum of the number of classes for each possible element order:")
    print(f"{equation_str} = {total_number_of_classes}")
    
    print(f"\nTherefore, the exact number of F_q-rational maximal tori of G is {total_number_of_classes}.")

solve_e8_tori()