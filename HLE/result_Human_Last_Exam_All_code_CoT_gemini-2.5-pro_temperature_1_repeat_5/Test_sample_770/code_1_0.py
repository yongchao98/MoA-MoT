import math

def calculate_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of C^3/I.

    According to the 3D McKay Correspondence, this rank is equal to the number
    of conjugacy classes of the icosahedral group I with age 1.
    """
    # The icosahedral group I is a subgroup of SO(3). Its conjugacy classes
    # are characterized by their rotation angles. There are 5 classes:
    # 1. Identity
    # 2. Rotation by pi (180 degrees) - axis through midpoints of opposite edges.
    # 3. Rotation by 2*pi/3 (120 degrees) - axis through centers of opposite faces.
    # 4. Rotation by 2*pi/5 (72 degrees) - axis through opposite vertices.
    # 5. Rotation by 4*pi/5 (144 degrees) - axis through opposite vertices.

    # We represent each class by its name and rotation angle phi.
    conjugacy_classes = {
        "Identity": 0,
        "Order 2 rotations": math.pi,
        "Order 3 rotations": 2 * math.pi / 3,
        "Order 5 rotations (type 1)": 2 * math.pi / 5,
        "Order 5 rotations (type 2)": 4 * math.pi / 5,
    }

    print("Calculating the 'age' for each conjugacy class of the icosahedral group:")
    print("-" * 70)

    age_one_classes_count = 0

    for name, phi in conjugacy_classes.items():
        print(f"Class: {name}")
        
        if phi == 0:
            # Identity element. Eigenvalues are (1, 1, 1).
            # Exponents (theta_j) are (0, 0, 0).
            age = 0
            print(f"  Rotation angle = 0 rad")
            print(f"  Eigenvalue exponents (theta_j) = (0.0, 0.0, 0.0)")
            print(f"  Age = 0.0 + 0.0 + 0.0 = {age}")
        else:
            # Non-identity rotation. Eigenvalues are (1, exp(i*phi), exp(-i*phi)).
            # We write them as exp(2*pi*i*theta_j) with 0 <= theta_j < 1.
            # The exponents theta_j are (0, phi/(2*pi), 1 - phi/(2*pi)).
            theta1 = 0
            theta2 = phi / (2 * math.pi)
            theta3 = 1 - theta2
            age = theta1 + theta2 + theta3
            
            print(f"  Rotation angle = {phi/math.pi:.2f}*pi rad")
            print(f"  Eigenvalue exponents (theta_j) = ({theta1:.2f}, {theta2:.2f}, {theta3:.2f})")
            # The final equation is formatted to show the sum.
            print(f"  Age = {theta1:.2f} + {theta2:.2f} + {theta3:.2f} = {age:.2f}")

        if round(age) == 1:
            age_one_classes_count += 1
            print("  This class has age 1.\n")
        else:
            print("  This class does not have age 1.\n")

    print("-" * 70)
    print(f"The number of conjugacy classes with age 1 is: {age_one_classes_count}")
    print("This number corresponds to the rank of H^2_c(Y, Q).")
    
    # The final answer must be an integer.
    return int(age_one_classes_count)

# Execute the calculation and print the final result in the required format.
final_answer = calculate_rank()
# The final result is the count of classes with age 1.
# The calculation shows this is 4.
# So, we output the final number.
# print(f'<<<{final_answer}>>>')
# Let's just output the final answer directly as requested by the format.

print("\nFinal Answer:")
print(f"<<<{final_answer}>>>")
