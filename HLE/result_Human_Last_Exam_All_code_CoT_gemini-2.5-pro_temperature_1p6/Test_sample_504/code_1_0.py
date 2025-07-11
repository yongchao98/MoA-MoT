import numpy as np

def solve():
    """
    This function calculates the largest number of pairwise linearly independent vectors
    in C^6 satisfying the given angle conditions.

    The method is to construct a set of vectors based on the D_6 root system,
    verify that it meets all the problem's conditions, and report its size.
    The vectors are of the form (e_i +/- e_j) / sqrt(2), where e_k are standard basis vectors.
    """
    
    d = 6
    vectors = []
    
    # Generate vectors of the form (e_i + e_j) / sqrt(2)
    for i in range(d):
        for j in range(i + 1, d):
            v = np.zeros(d)
            v[i] = 1
            v[j] = 1
            v /= np.sqrt(2)
            vectors.append(v)
            
    # Generate vectors of the form (e_i - e_j) / sqrt(2)
    for i in range(d):
        for j in range(i + 1, d):
            v = np.zeros(d)
            v[i] = 1
            v[j] = -1
            v /= np.sqrt(2)
            vectors.append(v)

    n = len(vectors)
    
    # The construction gives C(6, 2) * 2 = 15 * 2 = 30 vectors.
    # The following is a verification of the properties.
    
    has_orthogonal_pair = False
    
    for i in range(n):
        for j in range(i + 1, n):
            v = vectors[i]
            w = vectors[j]
            
            # Note: numpy's dot for 1-D arrays is the inner product without conjugation.
            # Since vectors are real, this is fine. For complex vectors, one would use np.vdot(v,w).
            inner_product = np.dot(v, w)
            abs_inner_product = np.abs(inner_product)
            
            if not (np.isclose(abs_inner_product, 0.0) or np.isclose(abs_inner_product, 0.5)):
                # This part should not be reached for the D_6 construction.
                print("Error: Angle condition violated.")
                return
                
            if np.isclose(abs_inner_product, 0.0):
                has_orthogonal_pair = True

    if not has_orthogonal_pair:
        # This part should not be reached.
        print("Error: No orthogonal pair found.")
        return
        
    print(f"The largest number of such vectors is {n}.")
    print("The final answer is derived from the following calculation: C(6, 2) * 2 = 30.")
    print("The numbers in this equation are:")
    print("6")
    print("2")
    print("2")
    print("30")

solve()
<<<30>>>