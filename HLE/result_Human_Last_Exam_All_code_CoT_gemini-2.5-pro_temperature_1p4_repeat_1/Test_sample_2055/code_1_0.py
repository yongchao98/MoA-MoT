import math

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k), returning 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve():
    """
    Solves the problem by calculating the eigenvalues of the Seidel matrix
    and then finding the lcm of their absolute values.
    """
    v = 50
    k = 5

    # The graph vertices are k-subsets of a v-set.
    # Adjacency is for intersection size 3. In the Johnson scheme J(v,k),
    # this is the relation A_e where e = k - intersection_size = 5 - 3 = 2.
    # The eigenvalues of A_e are given by the Eberlein polynomials P_i(e).
    
    print("Step 1: Calculate the eigenvalues of the adjacency matrix A.")
    # Eigenvalues of A=A_2 are P_i(2) for i = 0, ..., k
    # P_i(2) = sum_{j=0 to 2} (-1)^j C(i,j) C(k-i, 2-j) C(v-k-i, 2-j)
    eigenvalues_A = []
    for i in range(k + 1):
        p_i_2 = (combinations(i, 0) * combinations(k - i, 2) * combinations(v - k - i, 2) -
                 combinations(i, 1) * combinations(k - i, 1) * combinations(v - k - i, 1) +
                 combinations(i, 2) * combinations(k - i, 0) * combinations(v - k - i, 0))
        eigenvalues_A.append(p_i_2)
        print(f"P_{i}(2) = {p_i_2}")

    # Step 2: Calculate the eigenvalues of the Seidel matrix S = J - I - 2A.
    print("\nStep 2: Calculate the eigenvalues of the Seidel matrix S.")
    eigenvalues_S = []
    
    # For i = 0, the eigenvalue of J is N = C(v,k), and for I is 1.
    n = combinations(v, k)
    theta_0 = n - 1 - 2 * eigenvalues_A[0]
    eigenvalues_S.append(theta_0)
    print(f"theta_0 = C(50,5) - 1 - 2 * P_0(2) = {n} - 1 - 2 * {eigenvalues_A[0]} = {theta_0}")

    # For i > 0, the eigenvalue of J is 0.
    for i in range(1, k + 1):
        theta_i = -1 - 2 * eigenvalues_A[i]
        eigenvalues_S.append(theta_i)
        print(f"theta_{i} = -1 - 2 * P_{i}(2) = -1 - 2 * {eigenvalues_A[i]} = {theta_i}")

    # Step 3: Compute the maximum order in the Smith group.
    # This is the lcm of the absolute values of the distinct eigenvalues of S.
    print("\nStep 3: Compute the maximum order, which is the lcm of the absolute values of the eigenvalues of S.")
    
    # Take absolute values for the lcm calculation
    abs_eigenvalues = [abs(e) for e in eigenvalues_S]
    
    # Compute the lcm of all eigenvalues. math.lcm can take multiple arguments since Python 3.9.
    max_order = math.lcm(*abs_eigenvalues)
    
    # Format the final equation string
    equation_str = "lcm(" + ", ".join(map(str, abs_eigenvalues)) + ")"
    print(f"The final equation is: {equation_str} = {max_order}")
    
    return max_order

# Execute the solution
final_answer = solve()
print(f"\nThe maximum order among all elements of the Smith group of S is {final_answer}.")
print(f"\n<<<{final_answer}>>>")
