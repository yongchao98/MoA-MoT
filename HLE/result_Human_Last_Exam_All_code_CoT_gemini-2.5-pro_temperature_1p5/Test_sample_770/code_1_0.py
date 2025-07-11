def calculate_cohomology_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution of C^3/G,
    where G is the icosahedral group (A_5).
    """

    # Step 1: State the theoretical foundation.
    print("The problem is to find the rank of H^2_c(Y, Q).")
    print("This rank is determined by the structure of the group G, the icosahedral group, which is isomorphic to A_5.")
    print("\nAccording to the McKay correspondence for 3D singularities and Poincaré duality, this rank is equal to:")
    print("rank(H^2_c(Y, Q)) = dim(H_4(Y, Q)) = number of exceptional divisors in the resolution Y.")
    print("This number, in turn, is equal to the number of non-trivial conjugacy classes of the group G = A_5.")
    
    # Step 2: Define and count the conjugacy classes of A_5.
    # The order of A_5 is 60.
    conjugacy_classes_A5 = {
        "Identity (1)": 1,
        "3-cycles (e.g., (123))": 20,
        "Double transpositions (e.g., (12)(34))": 15,
        "5-cycles Class 1 (e.g., (12345))": 12,
        "5-cycles Class 2 (e.g., (12354))": 12
    }
    
    num_total_classes = len(conjugacy_classes_A5)
    
    print(f"\nThe group A_5 has {num_total_classes} conjugacy classes in total:")
    for class_name, size in conjugacy_classes_A5.items():
        print(f"- {class_name}: contains {size} elements.")
        
    # Step 3: Calculate the number of non-trivial classes.
    num_trivial_classes = 1 # The class of the identity element
    num_nontrivial_classes = num_total_classes - num_trivial_classes

    print("\nThe number of non-trivial conjugacy classes is the total number of classes minus the trivial (identity) class.")
    print("\nFinal Calculation:")
    print(f"{num_total_classes} - {num_trivial_classes} = {num_nontrivial_classes}")

    print(f"\nThus, the rank of H^2_c(Y, Q) is {num_nontrivial_classes}.")

calculate_cohomology_rank()
