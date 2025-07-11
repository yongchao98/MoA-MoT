import numpy as np

def solve():
    """
    Solves the problem based on the derived properties of the matrices.
    """
    
    # Step 1-3: Determine the function f(n).
    # Based on our analysis, the eigenvalues of W_n are all 1.
    # Therefore, f(n), the sum of the absolute cubes of the eigenvalues, is simply n.
    # f(n) = sum(|1|^3 for _ in range(n)) = n.
    
    # Step 4: Find the smallest n where f(n) > 10.
    # We need to find the smallest integer n such that n > 10.
    n = 10
    while True:
        # In our case f(n) = n, so we loop until n > 10.
        f_n = n 
        if f_n > 10:
            break
        n += 1
        
    print(f"The smallest integer n where f(n) > 10 is: {n}")

    # Step 5: Determine W_n and its infinity norm.
    # For this n, W_n is the Jordan block J_n(1).
    # We construct this matrix to calculate its infinity norm.
    # It has 1s on the diagonal and 1s on the superdiagonal.
    W_n = np.eye(n) + np.diag(np.ones(n - 1), 1)
    
    # The infinity norm is the maximum absolute row sum.
    norm_W_n = np.linalg.norm(W_n, np.inf)
    
    print(f"The infinity norm of W_{n} is: {norm_W_n}")

    # Step 6: Calculate the final result.
    result = n * norm_W_n
    
    print(f"\nThe final equation is: n * ||W_n||_inf = {n} * {norm_W_n} = {result}")
    
    # The problem also asks to return the answer in a specific format
    # This part is for the final answer submission and is not printed to the console.
    return result

# Execute the solution
if __name__ == "__main__":
    solve()
