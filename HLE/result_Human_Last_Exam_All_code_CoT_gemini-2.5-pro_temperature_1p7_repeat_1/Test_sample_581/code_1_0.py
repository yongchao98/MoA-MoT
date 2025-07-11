def solve_cap_set_bound():
    """
    Calculates lower bounds for the cap set problem in dimension 8
    using the product construction and states the best-known lower bound.
    """
    print("The cap set problem involves finding the maximum number of points in an n-dimensional vector space over the field of 3 elements, such that no three points are collinear.")
    print("We can find lower bounds for the size of a cap set in dimension n by using the product construction.")
    print("This means the size of a cap set in dimension a+b is at least the product of the sizes of cap sets in dimensions a and b.")
    print("\nLet's calculate some lower bounds for dimension 8 using known results from smaller dimensions:")

    # Known maximal sizes (or best lower bounds) for cap sets in dimension n
    cap_set_sizes = {
        2: 4,
        3: 9,
        4: 20,
        5: 45,
        6: 112
    }
    
    # Construction from dimension 4 + 4
    dim_a1, dim_a2 = 4, 4
    size_a1, size_a2 = cap_set_sizes[dim_a1], cap_set_sizes[dim_a2]
    bound_a = size_a1 * size_a2
    print(f"\nUsing dimensions {dim_a1} + {dim_a2}:")
    print(f"The equation is: {size_a1} * {size_a2} = {bound_a}")

    # Construction from dimension 5 + 3
    dim_b1, dim_b2 = 5, 3
    size_b1, size_b2 = cap_set_sizes[dim_b1], cap_set_sizes[dim_b2]
    bound_b = size_b1 * size_b2
    print(f"\nUsing dimensions {dim_b1} + {dim_b2}:")
    print(f"The equation is: {size_b1} * {size_b2} = {bound_b}")
    
    # Construction from dimension 6 + 2
    dim_c1, dim_c2 = 6, 2
    size_c1, size_c2 = cap_set_sizes[dim_c1], cap_set_sizes[dim_c2]
    bound_c = size_c1 * size_c2
    print(f"\nUsing dimensions {dim_c1} + {dim_c2}:")
    print(f"The equation is: {size_c1} * {size_c2} = {bound_c}")
    
    print("\nThe product constructions show a lower bound of at least 448.")
    print("However, more advanced constructions exist.")
    
    best_known_lower_bound = 496
    print(f"\nThe best known lower bound, from a 2008 construction by Yves Edel, is {best_known_lower_bound}.")

solve_cap_set_bound()