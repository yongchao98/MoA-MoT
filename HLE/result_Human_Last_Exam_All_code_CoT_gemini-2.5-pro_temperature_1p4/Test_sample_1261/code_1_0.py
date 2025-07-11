import numpy as np

def solve_binary_system_mod_q(n, k, m_factor=1.0):
    """
    Sets up a small instance of the problem Ax=0 (mod q) and finds a solution.
    
    The condition is m = Omega(n^k). For this demonstration, we'll use m = ceil(m_factor * n^k).
    Since the algorithm is theoretical, the code will find a solution via brute-force search, 
    which is feasible only for small m.
    """
    if n**k > 20:
        print("Parameters n and k are too large for a demonstration via brute-force search.")
        return

    q = 2**k
    # We set m to satisfy the Omega(n^k) condition for our small case.
    m = int(np.ceil(m_factor * (n**k)))
    if m <= n * k:
        print(f"Warning: m={m} may not be large enough compared to n*k={n*k}. Increasing m.")
        m = n * k + 1

    print(f"Goal: Find a non-zero vector x in {{0,1}}^{m} such that Ax = 0 (mod q)")
    print(f"Parameters: n={n}, k={k}, q=2^{k}={q}, m={m}\n")

    # 1. Generate a random matrix A from Z_q^{n x m}
    A = np.random.randint(0, q, size=(n, m))
    print("Randomly generated matrix A (mod q):")
    print(A)
    print("-" * 20)

    # 2. Search for a non-zero binary vector x
    # We iterate through all possible non-zero binary vectors x of length m.
    # The number of such vectors is 2^m - 1.
    solution_found = False
    # Start from 1, since we want a non-zero vector x.
    for i in range(1, 2**m):
        # Construct the binary vector x
        binary_representation = bin(i)[2:].zfill(m)
        x = np.array([int(bit) for bit in binary_representation])

        # Calculate Ax (mod q)
        # Using integer arithmetic then taking the modulus
        result_vector = (A @ x) % q

        # Check if Ax is the zero vector
        if np.all(result_vector == 0):
            print("Solution found!")
            print(f"Non-zero binary vector x:")
            print(x)
            print("\nVerification: Ax (mod q):")
            print(f"{A} @ {x} = {A @ x}")
            print(f"({A @ x}) mod {q} = {(A @ x) % q}")
            
            print("\nFinal Equation:")
            for row_idx in range(n):
                equation_parts = []
                for col_idx in range(m):
                    if x[col_idx] == 1:
                        equation_parts.append(f"{A[row_idx, col_idx]}*{x[col_idx]}")
                print(f"(Row {row_idx+1})  {' + '.join(equation_parts)} = {np.dot(A[row_idx, :], x)} \u2261 0 (mod {q})")

            solution_found = True
            break
            
    if not solution_found:
        print("No solution found for this random matrix A.")
        print("Note: A solution is guaranteed to exist for large enough n, but might not for small toy examples.")


# Let's choose small parameters that satisfy the problem conditions.
# Let k=2 (since k>1).
# Let n=2.
# Then m must be Omega(n^k) = Omega(2^2) = Omega(4). Let's choose m=5.
solve_binary_system_mod_q(n=2, k=2)
