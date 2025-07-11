import math

def solve_cohomology_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of the
    icosahedral quotient singularity C^3/G.

    The rank is determined by counting the number of "junior" (age 1)
    conjugacy classes of the icosahedral group G = A_5.
    """

    # The 5 conjugacy classes of the icosahedral group (A_5) are represented
    # by their corresponding rotation angles in the standard 3D representation.
    classes = [
        {"name": "Identity", "order": 1, "angle_in_rad": 0},
        {"name": "Rotation by pi", "order": 2, "angle_in_rad": math.pi},
        {"name": "Rotation by 2*pi/3", "order": 3, "angle_in_rad": 2 * math.pi / 3},
        {"name": "Rotation by 2*pi/5", "order": 5, "angle_in_rad": 2 * math.pi / 5},
        {"name": "Rotation by 4*pi/5", "order": 5, "angle_in_rad": 4 * math.pi / 5},
    ]

    print("Step-by-step calculation of the 'age' for each conjugacy class:")
    print("=" * 75)

    junior_class_count = 0

    for c in classes:
        name = c["name"]
        theta = c["angle_in_rad"]

        # For a rotation by angle theta, the eigenvalues are 1, exp(i*theta), exp(-i*theta).
        # We write them as exp(2*pi*i*alpha) with alpha in [0, 1).
        # The age is the sum of these alpha values.

        # For eigenvalue 1 = exp(2*pi*i*0):
        alpha1 = 0
        # For eigenvalue exp(i*theta) = exp(2*pi*i * theta/(2*pi)):
        alpha2 = theta / (2 * math.pi)
        # For eigenvalue exp(-i*theta) = exp(2*pi*i * (1 - theta/(2*pi))):
        # Note: if theta is 0, this alpha is also 0.
        alpha3 = (2 * math.pi - theta) / (2 * math.pi) if theta != 0 else 0

        # The age is the sum of the alpha values.
        age = alpha1 + alpha2 + alpha3

        print(f"Class: {name}")
        print(f"  - Rotation angle (theta): {theta/math.pi:.2f} * pi radians")
        print(f"  - Eigenvalue exponents (alpha_j): {alpha1:.3f}, {alpha2:.3f}, {alpha3:.3f}")
        # Final equation for the age calculation
        print(f"  - Age = {alpha1:.3f} + {alpha2:.3f} + {alpha3:.3f} = {age:.3f}")

        if round(age) == 1:
            junior_class_count += 1
            print("  - Result: This is a junior class (age = 1).")
        elif round(age) == 0:
            print("  - Result: This is the identity class (age = 0).")
        else:
            print(f"  - Result: This is not a junior class (age = {age:.3f}).")
        print("-" * 75)

    print(f"\nSummary: The number of junior (age 1) conjugacy classes is {junior_class_count}.")
    print("According to the McKay correspondence, this is the rank of H^2_c(Y, Q).")

solve_cohomology_rank()
<<<4>>>