import math

def solve_cohomology_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of C^3/A_5.
    
    The calculation follows these steps:
    1. The rank of H^2_c(Y, Q) is equal to b_4(Y) by Poincaré Duality.
    2. b_4(Y) is equal to b_2(Y), also by Poincaré Duality for Q-homology manifolds.
    3. b_2(Y) is the number of age-1 conjugacy classes of the group A_5.
    4. For any g in A_5 (subset SO(3)), age(g) = 1 if g is not the identity.
    5. The problem reduces to counting the non-identity conjugacy classes of A_5.
    """

    # A5 has 60 elements. Its conjugacy classes are well-known.
    # We list them here with their geometric interpretation in SO(3).
    A5_conjugacy_classes = [
        {'name': 'Identity', 'order': 1, 'size': 1, 'rotation_angle_rad': 0},
        {'name': 'Rotation by pi', 'order': 2, 'size': 15, 'rotation_angle_rad': math.pi},
        {'name': 'Rotation by 2pi/3', 'order': 3, 'size': 20, 'rotation_angle_rad': 2 * math.pi / 3},
        {'name': 'Rotation by 2pi/5', 'order': 5, 'size': 12, 'rotation_angle_rad': 2 * math.pi / 5},
        {'name': 'Rotation by 4pi/5', 'order': 5, 'size': 12, 'rotation_angle_rad': 4 * math.pi / 5}
    ]

    num_age_1_classes = 0
    
    print("Calculating the age for each conjugacy class of A5 in SO(3):")
    for C in A5_conjugacy_classes:
        angle = C['rotation_angle_rad']
        if angle == 0:
            # Identity element. Eigenvalues are (1, 1, 1).
            # The exponents theta_j are (0, 0, 0).
            age = 0
        else:
            # Non-identity element. Eigenvalues are (1, exp(i*theta), exp(-i*theta)).
            # The exponents theta_j are (0, theta/(2*pi), 1 - theta/(2*pi)).
            theta1 = 0
            theta2 = angle / (2 * math.pi)
            theta3 = 1 - theta2
            age = theta1 + theta2 + theta3
        
        print(f"- Class: {C['name']}, Order: {C['order']}, Size: {C['size']}")
        print(f"  Rotation angle: {angle/math.pi:.2f}*pi rad. Age = {age:.1f}")
        
        if age == 1:
            num_age_1_classes += 1

    b2_Y = num_age_1_classes
    b4_Y = b2_Y  # From Poincaré Duality b_4 = b_2
    rank_H2c = b4_Y

    print("\n--- Final Calculation ---")
    print("The rank of the second Betti number, b2(Y), is the number of age-1 conjugacy classes.")
    print(f"b2(Y) = {b2_Y}")
    print("\nBy Poincaré Duality for rational homology manifolds, b4(Y) = b2(Y).")
    print(f"b4(Y) = {b2_Y} = {b4_Y}")
    print("\nBy Poincaré Duality for non-compact manifolds, rank(H^2_c(Y, Q)) = b4(Y).")
    print(f"rank(H^2_c(Y, Q)) = b4(Y) = {rank_H2c}")

    return rank_H2c

if __name__ == '__main__':
    final_answer = solve_cohomology_rank()
    print(f"\nFinal Answer: {final_answer}")
    print(f'<<<4>>>')