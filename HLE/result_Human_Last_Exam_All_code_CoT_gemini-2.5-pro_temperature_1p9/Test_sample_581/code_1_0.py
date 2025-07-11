def find_cap_set_lower_bound():
    """
    This function calculates and explains the best-known lower bound for a cap set in dimension 8.
    """
    
    # The size of the largest cap set in dimension n over F_3 is denoted by r_3(n).
    # These values are known exactly for n <= 6.
    known_r3_values = {
        1: 2,
        2: 4,
        3: 9,
        4: 20,
        5: 45,
        6: 112,
    }

    target_dimension = 8
    
    print(f"The task is to find the best known lower bound for the size of a cap set in dimension {target_dimension}.")
    print("A common method to find a lower bound is the 'product construction'.")
    print(f"This method uses the sizes of cap sets in smaller dimensions that sum up to {target_dimension}.\n")

    best_product_bound = 0
    
    # Calculate the bound for the partition 2 + 6
    dim1, dim2 = 2, 6
    val1, val2 = known_r3_values[dim1], known_r3_values[dim2]
    product_bound_2_6 = val1 * val2
    print(f"For the partition {dim1} + {dim2}, the lower bound is r_3({dim1}) * r_3({dim2}):")
    print(f"{val1} * {val2} = {product_bound_2_6}")
    if product_bound_2_6 > best_product_bound:
        best_product_bound = product_bound_2_6
        
    # Calculate the bound for the partition 3 + 5
    dim1, dim2 = 3, 5
    val1, val2 = known_r3_values[dim1], known_r3_values[dim2]
    product_bound_3_5 = val1 * val2
    print(f"\nFor the partition {dim1} + {dim2}, the lower bound is r_3({dim1}) * r_3({dim2}):")
    print(f"{val1} * {val2} = {product_bound_3_5}")
    if product_bound_3_5 > best_product_bound:
        best_product_bound = product_bound_3_5

    # Calculate the bound for the partition 4 + 4
    dim1, dim2 = 4, 4
    val1, val2 = known_r3_values[dim1], known_r3_values[dim2]
    product_bound_4_4 = val1 * val2
    print(f"\nFor the partition {dim1} + {dim2}, the lower bound is r_3({dim1}) * r_3({dim2}):")
    print(f"{val1} * {val2} = {product_bound_4_4}")
    if product_bound_4_4 > best_product_bound:
        best_product_bound = product_bound_4_4

    print(f"\nThe best lower bound from these product constructions is {best_product_bound}.\n")
    
    print("However, more advanced algebraic constructions can produce larger sets.")
    print("In 2000, Yves Edel and Jurgen Bierbrauer developed a construction for dimension 8.")
    
    best_known_lower_bound = 496
    print(f"Their construction yields a cap set of size {best_known_lower_bound}.")
    print("\nThis is the best-known lower bound for the size of a cap set in dimension 8 to date.")

find_cap_set_lower_bound()