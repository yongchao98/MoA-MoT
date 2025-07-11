import math

def solve_h3_regular_elements():
    """
    Calculates the number of elements of the reflection group H3
    that have a regular eigenvector with an eigenvalue of order 10.
    """

    # Step 1: Define properties of the H3 reflection group
    group_name = "H3"
    group_order = 120  # |H3|, the order of the symmetry group of the icosahedron
    coxeter_number = 10  # h, the Coxeter number for H3
    exponents = [1, 5, 9]  # The exponents m_i for H3

    print(f"Properties of the reflection group of type {group_name}:")
    print(f"  - Group Order |W| = {group_order}")
    print(f"  - Coxeter Number h = {coxeter_number}")
    print(f"  - Exponents m_i = {exponents}\n")

    # Step 2: The question asks for elements with a regular eigenvector.
    # By definition, these are the "regular elements" of the group.
    # Regular elements are precisely those in the conjugacy class of a Coxeter element.
    print("An element has a regular eigenvector if and only if it is a 'regular element'.")
    print("The number of regular elements is the size of the conjugacy class of a Coxeter element.\n")

    # Step 3: Calculate the number of regular elements.
    # The size of this class is |W| / h.
    num_regular_elements = group_order // coxeter_number
    print("Calculating the number of regular elements:")
    print(f"  Number = |W| / h")
    print(f"  Number = {group_order} / {coxeter_number} = {num_regular_elements}\n")

    # Step 4: Verify these elements have an eigenvalue of order 10 with a regular eigenvector.
    print("Verifying the eigenvalue conditions for these regular elements:")
    print("The eigenvalues of a regular element are of the form exp(2*pi*i * m/h).")
    print("The eigenvector for an exponent 'm' is regular if and only if gcd(m, h) = 1.")
    print("The order of the corresponding eigenvalue is h / gcd(m, h), which is h=10 if and only if gcd(m, h) = 1.\n")

    print("Checking this condition for the exponents of H3:")
    found_satisfying_eigenvalue = False
    for m in exponents:
        g = math.gcd(m, coxeter_number)
        order = coxeter_number // g
        is_regular = (g == 1)
        if order == coxeter_number:
            print(f"  - For exponent m = {m}:")
            print(f"    gcd({m}, {coxeter_number}) = {g}. The eigenvector is regular.")
            print(f"    The eigenvalue order is {coxeter_number} / {g} = {order}.")
            found_satisfying_eigenvalue = True
        else:
            print(f"  - For exponent m = {m}:")
            print(f"    gcd({m}, {coxeter_number}) = {g}. The eigenvector is NOT regular.")
            print(f"    The eigenvalue order is {coxeter_number} / {g} = {order}.")
    
    print("\nSince there are exponents satisfying the condition (m=1 and m=9), every regular element has")
    print("at least one regular eigenvector with an eigenvalue of order 10.")
    
    print("\nConclusion: The number of such elements is equal to the total number of regular elements.\n")
    print("Final Calculation:")
    print(f"Number of elements = {group_order} / {coxeter_number} = {num_regular_elements}")


solve_h3_regular_elements()