import math
import functools

def get_adj_eigenvalues(v, k, i):
    """
    Calculates the eigenvalues of the adjacency matrix A_i of the Johnson scheme J(v,k).
    Adjacency is defined by intersection size k-i.
    """
    eigenvalues = []
    for j in range(k + 1):
        lambda_j = 0
        for l in range(i + 1):
            try:
                term = ((-1)**l *
                        math.comb(j, l) *
                        math.comb(k - j, i - l) *
                        math.comb(v - k - j, i - l))
                lambda_j += term
            except ValueError:
                # This can happen if a combination is not defined (e.g., k < 0),
                # which should be treated as 0.
                pass
        eigenvalues.append(lambda_j)
    return eigenvalues

def get_seidel_eigenvalues(v, k, adj_eigenvalues):
    """
    Calculates the eigenvalues of the Seidel matrix S = J - I - 2A.
    """
    seidel_eigenvalues = []
    num_vertices = math.comb(v, k)
    
    # Eigenvalue for the all-ones vector
    lambda_0 = adj_eigenvalues[0]
    theta_0 = num_vertices - 1 - 2 * lambda_0
    seidel_eigenvalues.append(theta_0)
    
    # Other eigenvalues
    for i in range(1, len(adj_eigenvalues)):
        lambda_i = adj_eigenvalues[i]
        theta_i = -1 - 2 * lambda_i
        seidel_eigenvalues.append(theta_i)
        
    return seidel_eigenvalues

def gcd(a, b):
    """Computes the greatest common divisor of two integers."""
    while b:
        a, b = b, a % b
    return abs(a)

def lcm(a, b):
    """Computes the least common multiple of two integers."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // gcd(a, b)

def main():
    """
    Main function to solve the problem by calculating eigenvalues and their LCM.
    """
    # Parameters of the graph
    v = 50  # size of the base set
    k = 5   # size of the subsets (which are the vertices)
    
    # Adjacency is based on intersection size 3.
    # In the Johnson scheme J(v,k), the matrix A_i corresponds to intersection size k-i.
    # So, 3 = 5 - i gives i = 2.
    intersection_size_minus_k = 2
    
    # Step 1: Calculate eigenvalues of the adjacency matrix A
    adj_eigenvalues = get_adj_eigenvalues(v, k, intersection_size_minus_k)
    
    # Step 2: Calculate eigenvalues of the Seidel matrix S
    seidel_eigenvalues = get_seidel_eigenvalues(v, k, adj_eigenvalues)
    
    # Step 3: The max order is the LCM of the absolute values of the eigenvalues of S.
    abs_eigenvalues = [abs(val) for val in seidel_eigenvalues]
    
    # Calculate the final result by finding the LCM of all values.
    result = functools.reduce(lcm, abs_eigenvalues)

    print("The distinct eigenvalues of the Seidel matrix S are:")
    for val in seidel_eigenvalues:
        print(f"  {val}")
    print("\nThe maximum order of an element in the Smith group of S is the")
    print("least common multiple (LCM) of the absolute values of these eigenvalues.")
    
    # Format the final equation string
    equation_str = ", ".join(map(str, abs_eigenvalues))
    print(f"\nWe need to compute:\nlcm({equation_str})")
    
    print(f"\nThe result is: {result}")

if __name__ == '__main__':
    main()