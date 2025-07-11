def solve_artin_torsion():
    """
    Calculates the number of minimal length positive torsion elements of a given
    order in the quotient of the E8 Artin group by its center.
    """
    # The problem asks for the number of torsion elements of order 10 in A/Z
    # for the Artin group A of type E8, which are positive words of minimal length.
    #
    # According to a theorem by Bessis, Digne, and Michel, the minimal length
    # for a positive torsion element of order 'd' is 'd' itself.
    # The number of such elements is given by the formula: d * N_d, where N_d
    # is the number of fundamental invariant degrees of the corresponding
    # Coxeter group W(E8) that are divisible by d.

    # Step 1: Define the degrees of the fundamental invariants for W(E8).
    # These are well-known mathematical constants associated with the E8 group.
    degrees_E8 = [2, 8, 12, 14, 18, 20, 24, 30]
    print("Step 1: The fundamental invariant degrees for the Coxeter group W(E8) are:")
    print(degrees_E8)
    print("-" * 30)

    # Step 2: Define the required order 'd'.
    d = 10
    print(f"Step 2: The required order of the torsion elements is d = {d}.")
    print("-" * 30)

    # Step 3: Find which of the degrees are divisible by d.
    divisible_degrees = [degree for degree in degrees_E8 if degree % d == 0]
    print(f"Step 3: The degrees from the list that are divisible by {d} are:")
    print(divisible_degrees)
    print("-" * 30)

    # Step 4: Count the number of these divisible degrees. This gives N_d.
    num_classes = len(divisible_degrees)
    print(f"Step 4: The number of such degrees, N_d, is {num_classes}.")
    print("-" * 30)

    # Step 5: Calculate the total number of elements using the formula d * N_d.
    num_elements = d * num_classes
    print("Step 5: The total number of elements is calculated using the formula d * N_d.")
    print(f"The final equation is: {d} * {num_classes} = {num_elements}")
    print("-" * 30)
    
    print(f"Therefore, there are {num_elements} torsion elements of order 10 with minimal positive word length.")

solve_artin_torsion()