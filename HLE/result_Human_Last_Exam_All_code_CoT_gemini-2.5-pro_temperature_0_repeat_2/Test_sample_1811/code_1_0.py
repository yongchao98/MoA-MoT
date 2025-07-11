def solve_vector_field_zeros():
    """
    Explains and prints the solution to the problem of finding the minimum
    number of zeros of a vector field on a compact manifold with boundary.
    """
    print("The problem asks for the least number of zeros a vector field can have on a compact manifold M with a non-empty boundary.")
    print("The answer is determined entirely by the Euler characteristic of the manifold, chi(M).")
    print("\n--- The Key Theorem ---")
    print("A fundamental result in differential topology states that a compact manifold M admits a vector field with no zeros if and only if its Euler characteristic, chi(M), is zero.")
    
    print("\n--- Analysis of Cases ---")
    
    print("\nCase 1: chi(M) = 0")
    print("According to the theorem, if chi(M) = 0, a vector field with 0 zeros can be constructed.")
    print("Therefore, the minimum number of zeros is 0.")
    
    print("\nCase 2: chi(M) != 0")
    print("If chi(M) is not 0, the theorem implies that any vector field must have at least one zero. So the minimum is at least 1.")
    print("Furthermore, it is always possible to construct a vector field with exactly one zero. This is done by creating a vector field where the sum of the indices of its zeros equals chi(M), and then mathematically 'merging' all the zeros into a single point.")
    print("Since a vector field with one zero can be constructed, the minimum number of zeros is 1.")
    
    print("\n--- Conclusion ---")
    print("The least number of zeros depends only on chi(M) and not on the Euler characteristic of its boundary, chi(∂M).")
    print("The final answer is given by the following rule:")
    print("\n*****************************************")
    print("* Least number of zeros = 0, if chi(M) = 0  *")
    print("* Least number of zeros = 1, if chi(M) != 0 *")
    print("*****************************************")

if __name__ == '__main__':
    solve_vector_field_zeros()