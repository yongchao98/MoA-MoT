import math
from fractions import Fraction

def solve_cohomology_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of C^3/G,
    where G is the icosahedral group (A_5).

    The methodology is based on the 3-dimensional McKay correspondence.
    """
    print("Step 1: Understand the theoretical background.")
    print("The rank of H^2_c(Y, Q) is equal to the Betti number b_2(Y).")
    print("The McKay correspondence for 3-folds states that b_2(Y) equals the")
    print("number of conjugacy classes of the group G with age equal to 1.\n")

    print("Step 2: Define the group G and its conjugacy classes.")
    print("G is the icosahedral group, isomorphic to A_5. Its action on C^3")
    print("is its standard 3D representation in SO(3).")
    
    # Conjugacy classes are defined by their order and corresponding rotation angle theta.
    # A rotation by theta has eigenvalues 1, exp(i*theta), exp(-i*theta).
    conjugacy_classes = [
        ("Identity", 1, 0),
        ("Rotations of order 3", 3, 2 * math.pi / 3),
        ("Rotations of order 2 (180-deg)", 2, math.pi),
        ("Rotations of order 5 (72-deg)", 5, 2 * math.pi / 5),
        ("Rotations of order 5 (144-deg)", 5, 4 * math.pi / 5)
    ]
    print(f"There are {len(conjugacy_classes)} conjugacy classes in A_5.\n")
    
    print("Step 3: Calculate the 'age' for each conjugacy class.")
    print("Age is sum(phi_j) where eigenvalues are exp(2*pi*i*phi_j) with phi_j in [0, 1).")
    print("-" * 60)

    age_1_classes_count = 0
    equation_terms = []

    for name, order, theta in conjugacy_classes:
        print(f"Analyzing class: {name}")
        
        # Calculate the exponents phi_j for the eigenvalues.
        if theta == 0:
            # For the identity element
            phi1 = Fraction(0)
            phi2 = Fraction(0)
            phi3 = Fraction(0)
        else:
            # For a rotation by theta
            phi1 = Fraction(0)  # For eigenvalue 1
            phi2 = Fraction(theta / (2 * math.pi)).limit_denominator() # For exp(i*theta)
            # For exp(-i*theta), which corresponds to exp(i*(2*pi-theta))
            phi3 = Fraction((2 * math.pi - theta) / (2 * math.pi)).limit_denominator()

        # The age is the sum of the exponents.
        age = phi1 + phi2 + phi3

        print(f"  Eigenvalue exponents (phi_j): {phi1}, {phi2}, {phi3}")
        print(f"  Age = {phi1} + {phi2} + {phi3} = {age}")

        if age == 1:
            age_1_classes_count += 1
            equation_terms.append("1")
            print("  This class contributes to the rank.")
        else:
            print(f"  This class has age {age} and does not contribute.")
        print("-" * 60)

    print("Step 4: Sum the contributions to find the final rank.")
    final_equation = " + ".join(equation_terms)
    print(f"The rank of H^2_c(Y, Q) is the total count of classes with age 1.")
    print(f"Rank = {final_equation} = {age_1_classes_count}")

solve_cohomology_rank()