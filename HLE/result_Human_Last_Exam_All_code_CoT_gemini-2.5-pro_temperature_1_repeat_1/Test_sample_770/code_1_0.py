import math

def calculate_age_and_rank():
    """
    Calculates the rank of H^2_c(Y, Q) by finding the number of
    age-1 conjugacy classes of the icosahedral group A_5.
    """
    # The conjugacy classes of A_5 correspond to types of rotations of an icosahedron.
    # We represent them by their order and the rotation angle.
    # Note: 5-cycles in S_5 split into two conjugacy classes in A_5,
    # corresponding to rotations by 2*pi/5 and 4*pi/5.
    
    conjugacy_classes = {
        "Identity (order 1)": 0,
        "Rotation by 2*pi/3 (order 3)": (2 * math.pi) / 3,
        "Rotation by pi (order 2)": math.pi,
        "Rotation by 2*pi/5 (order 5, class 1)": (2 * math.pi) / 5,
        "Rotation by 4*pi/5 (order 5, class 2)": (4 * math.pi) / 5,
    }

    print("Calculating the age for each conjugacy class of the icosahedral group (A_5) in its 3D representation.")
    print("-" * 80)
    
    age_one_classes = 0

    for name, angle in conjugacy_classes.items():
        # For a rotation in SO(3) by angle alpha, the eigenvalues are 1, e^(i*alpha), e^(-i*alpha).
        # The eigenvalues lambda_j are written as exp(2*pi*i*theta_j) with 0 <= theta_j < 1.
        
        if angle == 0:
            # Identity element
            theta1, theta2, theta3 = 0.0, 0.0, 0.0
        else:
            # For a rotation by `angle`, the theta values are 0, angle/(2*pi), (2*pi-angle)/(2*pi)
            # which simplifies to angle/(2*pi) and 1 - angle/(2*pi)
            theta1 = 0.0
            theta2 = angle / (2 * math.pi)
            theta3 = 1.0 - theta2
            
        age = theta1 + theta2 + theta3
        
        print(f"Class: {name}")
        # Using f-strings to format the output nicely to two decimal places
        print(f"  Rotation angle (radians): {angle:.2f}")
        print(f"  Theta values: {theta1:.2f}, {theta2:.2f}, {theta3:.2f}")
        print(f"  Age = {theta1:.2f} + {theta2:.2f} + {theta3:.2f} = {age:.2f}")

        if math.isclose(age, 1.0):
            age_one_classes += 1
            print("  This is an age-1 class.\n")
        else:
            print("  This is not an age-1 class.\n")

    print("-" * 80)
    print(f"The number of conjugacy classes with age 1 is {age_one_classes}.")
    print("This number corresponds to the number of exceptional divisors in the resolution Y,")
    print("which is the rank of H_2(Y, Q).")
    print("By Poincare duality, rank(H^2_c(Y, Q)) = rank(H^4(Y, Q)) = rank(H_2(Y, Q)).")
    print(f"Therefore, the rank of H^2_c(Y, Q) is {age_one_classes}.")

calculate_age_and_rank()
<<<4>>>