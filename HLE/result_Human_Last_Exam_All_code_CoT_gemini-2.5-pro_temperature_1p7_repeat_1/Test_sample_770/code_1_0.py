import numpy as np

def calculate_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of C^3/G,
    where G is the icosahedral group.
    """

    # The icosahedral group is isomorphic to A_5.
    # We list its conjugacy classes by their geometric interpretation (rotations),
    # their order, and the size of the class.
    # The identity class is non-junior (age=0). All non-trivial classes
    # in the standard 3D SO(3) representation have age=1.

    conjugacy_classes = [
        {"name": "Identity", "order": 1, "rotation_angle_rad": 0},
        {"name": "3-fold rotation", "order": 3, "rotation_angle_rad": 2 * np.pi / 3},
        {"name": "2-fold rotation", "order": 2, "rotation_angle_rad": np.pi},
        {"name": "5-fold rotation A", "order": 5, "rotation_angle_rad": 2 * np.pi / 5},
        {"name": "5-fold rotation B", "order": 5, "rotation_angle_rad": 4 * np.pi / 5},
    ]

    print("Step 1: Identify conjugacy classes of the icosahedral group (A_5) and calculate their 'age'.")
    print("The age of a group element g with eigenvalues e^(2*pi*i*theta_j) is sum(theta_j).")
    print("For a rotation by angle alpha in SO(3), the eigenvalues are 1, e^(i*alpha), e^(-i*alpha).")
    print("The corresponding phases theta_j in [0, 1) are 0, alpha/(2*pi), 1 - alpha/(2*pi).")
    print("The age is the sum of these phases.\n")

    junior_classes_count = 0
    for cc in conjugacy_classes:
        name = cc["name"]
        angle = cc["rotation_angle_rad"]

        if angle == 0:
            phases = [0, 0, 0]
        else:
            theta1 = 0
            theta2 = angle / (2 * np.pi)
            theta3 = 1 - theta2
            phases = [theta1, theta2, theta3]
        
        age = sum(phases)

        print(f"Class: {name}")
        print(f"  Rotation Angle: {angle/np.pi:.2f}*pi")
        print(f"  Phases (theta_j): [{phases[0]:.2f}, {phases[1]:.2f}, {phases[2]:.2f}]")
        print(f"  Age: {age:.0f}")

        # A class is junior if its age is 1.
        if age == 1:
            junior_classes_count += 1
            print("  This is a junior class.")
        else:
            print("  This is not a junior class (age=0).")
        print("-" * 20)
    
    b2_Y = junior_classes_count
    print(f"Step 2: According to the 3D McKay Correspondence, the second Betti number b2(Y) is the number of junior conjugacy classes.")
    print(f"Number of junior classes = {b2_Y}")
    print(f"So, b2(Y) = {b2_Y}")
    
    print("\nStep 3: Relate b2(Y) to rank H^2_c(Y, Q).")
    b4_Y = b2_Y
    rank_H2c = b4_Y
    
    print(f"For a crepant resolution Y of C^3/G, Y is a Calabi-Yau 3-fold.")
    print(f"By Poincaré Duality, rank H^2_c(Y, Q) = rank H^4(Y, Q) = b4(Y).")
    print(f"By Hodge duality on a Calabi-Yau 3-fold, b4(Y) = b2(Y).")
    print(f"Therefore, the final equation is:")
    print(f"rank H^2_c(Y, Q) = b4(Y) = b2(Y) = {junior_classes_count}")

    final_answer = rank_H2c
    return final_answer

if __name__ == '__main__':
    answer = calculate_rank()
    # The final answer will be printed inside the function as part of the explanation.
