def solve_h3_regular_elements():
    """
    Calculates the number of elements of the reflection group of type H3
    that have a regular eigenvector with an eigenvalue of order 10.
    """

    # The order of the reflection group of type H3.
    group_order_W = 120

    # The Coxeter number for the H3 group.
    coxeter_number_h = 10

    print("Step 1: Understanding the properties of the H3 group.")
    print("The problem asks for the number of elements with a regular eigenvector whose corresponding eigenvalue has order 10.")
    print("An element with a regular eigenvector is called a 'regular element'.")
    print("\\n")
    print("Step 2: Identifying the eigenvalues of regular elements in H3.")
    print("For any regular element in a Coxeter group, its eigenvalues are given by the formula: exp(2 * pi * i * m / h),")
    print("where 'h' is the Coxeter number and 'm' represents the exponents of the group.")
    print(f"For the H3 group, the Coxeter number h = {coxeter_number_h} and the exponents are 1, 5, and 9.")
    print("This means every regular element in H3 has eigenvalues with orders 10, 2, and 10.")
    print("Therefore, any regular element in H3 automatically satisfies the condition of having an eigenvalue of order 10.")
    print("\\n")
    print("Step 3: Calculating the number of regular elements.")
    print("The problem is now reduced to finding the total number of regular elements in H3.")
    print("The number of regular elements is given by the formula: |W| / h, where |W| is the order of the group.")
    print(f"For H3, |W| = {group_order_W} and h = {coxeter_number_h}.")
    print("\\n")

    # The number of such elements is the order of the group divided by the Coxeter number.
    num_elements = group_order_W // coxeter_number_h

    print("Step 4: The final calculation.")
    print(f"The number of elements is {group_order_W} / {coxeter_number_h} = {num_elements}.")


solve_h3_regular_elements()