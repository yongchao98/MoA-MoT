from fractions import Fraction

def solve_cohomology_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution of C^3/G,
    where G is the icosahedral group.

    The rank is equal to the number of age-1 conjugacy classes of G.
    """

    print("Starting the calculation for the rank of H^2_c(Y, Q).")
    print("This is equivalent to finding the number of conjugacy classes with age 1.\n")

    # The icosahedral group A_5 has 5 conjugacy classes.
    # Its action on C^3 is the standard one as a subgroup of SO(3).
    # A rotation in SO(3) by angle alpha has eigenvalues {1, exp(i*alpha), exp(-i*alpha)}.
    # The order of the element is n. The rotation angle alpha is a multiple of 2*pi/n.
    # For a rotation by alpha = 2*pi*k/n, the exponents r_j are {0, k/n, (n-k)/n}.
    # The age is their sum.

    conjugacy_classes = [
        {'order': 1, 'name': 'Identity', 'k': 0, 'size': 1},
        {'order': 2, 'name': 'Elements of order 2', 'k': 1, 'size': 15}, # Rotation by pi
        {'order': 3, 'name': 'Elements of order 3', 'k': 1, 'size': 20}, # Rotation by 2*pi/3
        {'order': 5, 'name': 'Elements of order 5 (type 1)', 'k': 1, 'size': 12}, # Rotation by 2*pi/5
        {'order': 5, 'name': 'Elements of order 5 (type 2)', 'k': 2, 'size': 12}, # Rotation by 4*pi/5
    ]

    age_1_class_count = 0

    for cc in conjugacy_classes:
        order = cc['order']
        name = cc['name']
        k = cc['k']

        print(f"--- Analyzing Class: {name} (Order {order}) ---")

        if order == 1:
            # The identity element has eigenvalues {1, 1, 1}
            r = [Fraction(0), Fraction(0), Fraction(0)]
        else:
            # Non-trivial rotations have eigenvalues {1, exp(2*pi*i*k/n), exp(-2*pi*i*k/n)}
            r = [Fraction(0), Fraction(k, order), Fraction(order - k, order)]

        age = sum(r)

        # Output the required equation with numbers
        num1, num2, num3 = float(r[0]), float(r[1]), float(r[2])
        sum_val = float(age)

        print(f"The eigenvalue exponents r_j are: [{r[0]}, {r[1]}, {r[2]}]")
        print(f"The equation for the age is: {num1:.3f} + {num2:.3f} + {num3:.3f} = {sum_val:.3f}")
        
        if age == 1:
            age_1_class_count += 1
            print("Result: This conjugacy class has age 1.\n")
        else:
            print(f"Result: This conjugacy class has age {age}.\n")

    print("==========================================================")
    print(f"The total number of conjugacy classes with age 1 is {age_1_class_count}.")
    print("This number corresponds to the rank of H^2_c(Y, Q).")
    print(f"Final Answer: {age_1_class_count}")
    return age_1_class_count

solve_cohomology_rank()
<<<4>>>