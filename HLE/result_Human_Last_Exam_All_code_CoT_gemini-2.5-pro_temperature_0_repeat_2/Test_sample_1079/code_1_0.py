def count_h3_elements():
    """
    Calculates the number of elements of the reflection group of type H3
    that have a regular eigenvector with a corresponding eigenvalue of order 10.
    """
    # Properties of the H3 reflection group
    # The order of the group H3 is 120.
    order_W = 120
    # The Coxeter number for H3 is 10.
    coxeter_number_h = 10

    print("The problem is to find the number of elements in H3 that are Coxeter elements.")
    print("This is because only Coxeter elements have regular eigenvectors, and for H3,")
    print("their corresponding eigenvalues are guaranteed to have order 10.")
    print("\nThe number of Coxeter elements is the order of the group divided by the Coxeter number.")
    
    # Calculate the number of Coxeter elements
    num_elements = order_W // coxeter_number_h
    
    print("\nFinal Calculation:")
    # The final output must show the numbers in the equation.
    print(f"{order_W} / {coxeter_number_h} = {num_elements}")

count_h3_elements()