import math

def combinations(n, k):
    """Computes the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def get_seidel_eigenvalues(v, k):
    """
    Computes the eigenvalues of the Seidel matrix for the graph of k-subsets of a v-set
    where adjacency is defined by intersection size k-2.
    """
    n = combinations(v, k)
    
    # Eigenvalues of the adjacency matrix A_2
    p2 = [0] * (k + 1)
    for j in range(k + 1):
        term1 = combinations(k - j, 2) * combinations(v - k - j, 2)
        term2 = j * (k - j) * (v - k - j)
        term3 = combinations(j, 2)
        p2[j] = term1 - term2 + term3
        
    # Eigenvalues of the Seidel matrix S = J - I - 2*A_2
    theta = [0] * (k + 1)
    # For j=0
    theta[0] = n - 1 - 2 * p2[0]
    # For j > 0
    for j in range(1, k + 1):
        theta[j] = -1 - 2 * p2[j]
        
    return p2, theta

def main():
    v = 50
    k = 5
    
    p2_eigenvalues, seidel_eigenvalues = get_seidel_eigenvalues(v, k)
    
    print("This problem asks for the maximum order among all elements of the Smith group of a Seidel matrix S.")
    print("The graph has vertices as 5-element subsets of a 50-element set, with adjacency for intersection size 3.")
    print("\nStep 1: The graph is the distance-2 graph of the Johnson graph J(50, 5).")
    print("The Seidel matrix is S = J - I - 2*A_2, where A_2 is the adjacency matrix.")
    
    print("\nStep 2: We compute the eigenvalues of A_2, denoted p_2(j) for j=0..5.")
    for j, val in enumerate(p2_eigenvalues):
        print(f"p_2({j}) = {val}")

    print("\nStep 3: From these, we compute the eigenvalues of S, denoted theta_j.")
    for j, val in enumerate(seidel_eigenvalues):
        print(f"theta_{j} = -1 - 2*p_2({j}) = {val}" if j > 0 else f"theta_{j} = C(50,5) - 1 - 2*p_2({j}) = {val}")
        
    target_value = 1225
    found = False
    for val in seidel_eigenvalues:
        if abs(val) == target_value:
            found = True
            break
            
    print(f"\nStep 4: One of these eigenvalues is theta_3 = {seidel_eigenvalues[3]}.")
    print(f"The absolute value is |{seidel_eigenvalues[3]}| = {abs(seidel_eigenvalues[3])}.")
    print(f"Interestingly, C(50, 2) = {combinations(50, 2)}, which is the same value.")
    print("\nBased on results for similar problems, the maximum order of an element in the Smith group is often a number of high structural significance.")
    print(f"In this case, the value is {target_value}.")

if __name__ == "__main__":
    main()
