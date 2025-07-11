import math

def combinations(n, k):
    """
    Calculates the number of combinations 'n choose k'.
    Returns 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def get_adjacency_eigenvalues(v, k):
    """
    Calculates the eigenvalues of the adjacency matrix A for the graph
    where vertices are k-subsets of a v-set and adjacency is intersection size k-2.
    """
    eigenvalues = []
    for j in range(k + 1):
        # Formula for eigenvalues of the 2nd associate of Johnson scheme J(v,k)
        term1 = combinations(k - j, 2) * combinations(v - k - j, 2)
        term2 = -j * (k - j) * (v - k - j)
        term3 = combinations(j, 2)
        theta_j = term1 + term2 + term3
        eigenvalues.append(theta_j)
    return eigenvalues

def get_seidel_eigenvalues(v, k, adj_eigenvalues):
    """
    Calculates the eigenvalues of the Seidel matrix S = J - I - 2A.
    """
    n = combinations(v, k)
    theta_0 = adj_eigenvalues[0]
    
    # Eigenvalue for the all-ones eigenvector
    sigma_0 = n - 1 - 2 * theta_0
    
    seidel_eigenvalues = [sigma_0]
    
    # Eigenvalues for vectors orthogonal to the all-ones vector
    for j in range(1, k + 1):
        theta_j = adj_eigenvalues[j]
        sigma_j = -1 - 2 * theta_j
        seidel_eigenvalues.append(sigma_j)
        
    return seidel_eigenvalues

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def lcm(a, b):
    """Computes the least common multiple of a and b."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // gcd(a, b)

def main():
    """
    Main function to solve the problem.
    """
    v = 50
    k = 5

    # Step 1 & 2: Get eigenvalues of the adjacency matrix A
    adj_eigenvalues = get_adjacency_eigenvalues(v, k)

    # Step 3: Get eigenvalues of the Seidel matrix S
    seidel_eigs = get_seidel_eigenvalues(v, k, adj_eigenvalues)
    
    distinct_eigs = sorted(list(set(seidel_eigs)))
    
    print("The distinct eigenvalues of the Seidel matrix S are:")
    for eig in distinct_eigs:
        print(eig)
        
    # Step 4 & 5: Calculate the LCM of the absolute values of the eigenvalues
    if not distinct_eigs:
        result = 0
    else:
        result = 1
        for eig in distinct_eigs:
            result = lcm(result, eig)

    # Output the final equation as requested
    equation_str = "lcm(" + ", ".join([f"|{eig}|" for eig in distinct_eigs]) + ")"
    print(f"\nThe maximum order is the LCM of the absolute values of these eigenvalues.")
    print(f"{equation_str} = {result}")


if __name__ == "__main__":
    main()
