def calculate_and_print_age(class_name, order, exponents, numerators):
    """
    Calculates the age for a given conjugacy class and prints the steps.
    
    Args:
        class_name (str): The name of the conjugacy class.
        order (int): The order of elements in the class.
        exponents (list of str): String representations of the fractional exponents.
        numerators (list of int): The numerators of the exponents for calculation.
    
    Returns:
        int: 1 if the age is 1, otherwise 0.
    """
    print(f"--- Analyzing {class_name} (order {order}) ---")
    
    # Calculate age
    age = sum(n / order for n in numerators)
    
    # Format the equation string
    equation = f"Age = {' + '.join(exponents)} = {age:.0f}"
    
    print(f"The eigenvalues correspond to fractional exponents {', '.join(exponents)}.")
    print(equation)
    
    if age == 1:
        print("This class has age 1.")
        print("-" * 40)
        return 1
    else:
        print(f"This class does not have age 1.")
        print("-" * 40)
        return 0

def solve_problem():
    """
    Solves the problem by calculating the number of age-1 conjugacy classes for A_5 in SL(3,C).
    """
    print("The rank of H^2_c(Y, Q) is the number of non-trivial conjugacy classes of G=A_5 with age 1.")
    print("We analyze the four non-trivial conjugacy classes of A_5.\n")
    
    age_1_class_count = 0
    
    # 1. Class of order 2 (e.g., permutation (12)(34))
    # Eigenvalues in the 3D irrep are {-1, -1, 1}.
    # As roots of unity, -1 = exp(2*pi*i * 1/2) and 1 = exp(2*pi*i * 0/2).
    # Exponents k/m are {1/2, 1/2, 0/2}.
    age_1_class_count += calculate_and_print_age(
        class_name="Class of order 2",
        order=2,
        exponents=["1/2", "1/2", "0/2"],
        numerators=[1, 1, 0]
    )

    # 2. Class of order 3 (e.g., permutation (123))
    # Eigenvalues are {1, w, w^2} where w = exp(2*pi*i/3).
    # Exponents k/m are {0/3, 1/3, 2/3}.
    age_1_class_count += calculate_and_print_age(
        class_name="Class of order 3",
        order=3,
        exponents=["0/3", "1/3", "2/3"],
        numerators=[0, 1, 2]
    )

    # 3. First class of order 5 (e.g., permutation (12345))
    # Eigenvalues are {1, z, z^4} where z = exp(2*pi*i/5).
    # Exponents k/m are {0/5, 1/5, 4/5}.
    age_1_class_count += calculate_and_print_age(
        class_name="First class of order 5",
        order=5,
        exponents=["0/5", "1/5", "4/5"],
        numerators=[0, 1, 4]
    )
    
    # 4. Second class of order 5 (e.g., permutation (13524))
    # Eigenvalues are {1, z^2, z^3} where z = exp(2*pi*i/5).
    # Exponents k/m are {0/5, 2/5, 3/5}.
    age_1_class_count += calculate_and_print_age(
        class_name="Second class of order 5",
        order=5,
        exponents=["0/5", "2/5", "3/5"],
        numerators=[0, 2, 3]
    )
    
    print(f"All four non-trivial conjugacy classes have age 1.")
    print(f"The total count is {age_1_class_count}.")
    print(f"\nThus, the rank of H^2_c(Y, Q) is {age_1_class_count}.")

solve_problem()