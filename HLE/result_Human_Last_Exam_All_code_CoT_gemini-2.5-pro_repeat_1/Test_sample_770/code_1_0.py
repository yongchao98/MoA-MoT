import math

def solve_cohomology_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of C^3/G,
    where G is the icosahedral group A_5.

    The rank of H^2_c(Y, Q) is equal to the second Betti number b_2(Y).
    By the 3D McKay Correspondence, b_2(Y) equals the number of conjugacy
    classes of G with age 1.

    For the standard action of G = A_5 in SO(3), any non-identity element g
    is a rotation. Its eigenvalues are {1, e^{i*theta}, e^{-i*theta}}.
    The age of g is the sum of the fractional parts of the eigenvalue phases,
    which is age(g) = 0 + theta/(2*pi) + (1 - theta/(2*pi)) = 1.

    Thus, the task is to count the non-identity conjugacy classes of A_5.
    """

    # The 5 conjugacy classes of the icosahedral group A_5,
    # characterized by element order and rotation angle for the standard 3D representation.
    conjugacy_classes = [
        {'name': '1A', 'order': 1, 'description': 'Identity', 'angle_rad': 0},
        {'name': '2A', 'order': 2, 'description': 'Rotation by pi', 'angle_rad': math.pi},
        {'name': '3A', 'order': 3, 'description': 'Rotation by 2pi/3', 'angle_rad': 2 * math.pi / 3},
        {'name': '5A', 'order': 5, 'description': 'Rotation by 2pi/5', 'angle_rad': 2 * math.pi / 5},
        {'name': '5B', 'order': 5, 'description': 'Rotation by 4pi/5', 'angle_rad': 4 * math.pi / 5},
    ]

    print("The rank of H^2_c(Y, Q) is the number of conjugacy classes with age 1.")
    print("Calculating the age for each conjugacy class of the icosahedral group A_5 acting on C^3.")
    print("-" * 60)

    junior_class_contributions = []

    for cc in conjugacy_classes:
        name = cc['name']
        order = cc['order']

        # The identity element has age 0.
        if order == 1:
            print(f"Class '{name}' (Identity): Age = 0. This class is not counted.")
            continue

        # For any non-identity element g in SO(3), which is a rotation by angle theta,
        # the eigenvalues are {1, exp(i*theta), exp(-i*theta)}.
        # The corresponding exponents a_j in [0,1) are {0, theta/(2*pi), 1-theta/(2*pi)}.
        # The age is the sum of these exponents, which is always 1.
        age = 1
        
        angle_deg = cc['angle_rad'] * 180 / math.pi
        print(f"Class '{name}': A rotation of order {order} by {angle_deg:.0f} degrees.")
        print(f"  > The age is calculated to be {age}.")
        junior_class_contributions.append("1")

    print("-" * 60)
    
    final_sum_str = " + ".join(junior_class_contributions)
    total_rank = len(junior_class_contributions)

    print("The total rank is the sum of contributions from each non-identity class:")
    print(f"Rank = {final_sum_str} = {total_rank}")

solve_cohomology_rank()
<<<4>>>