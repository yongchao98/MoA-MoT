import math

def calculate_curvature_cost(d, n):
    """
    Calculates the minimum curvature cost for the NGD update.

    Args:
        d (int): The dimension of the input/output space.
        n (int): The number of training samples, with n < d.
    """
    if n >= d:
        print("This analysis assumes n < d.")
        return

    print(f"Analysis for d = {d} and n = {n}:\n")

    # The NGD update cost can be reduced to O(d^2 * n) using the Woodbury identity.
    # The calculation is V = (1/alpha) * (G - (G @ X) @ inv(n*alpha*I + X.T @ X) @ X.T)
    # We analyze the cost of the most expensive operations in floating point operations (FLOPS).

    # Cost of X.T @ X: (n,d) @ (d,n) -> (n,n). Cost is d*n*n = d*n^2
    cost_xtx = d * n**2

    # Cost of inverting an (n,n) matrix: n^3
    cost_inv = n**3

    # Cost of G @ X: (d,d) @ (d,n) -> (d,n). Cost is d*d*n = d^2*n
    cost_gx = d**2 * n
    
    # Cost of (G@X) @ inv(M): (d,n) @ (n,n) -> (d,n). Cost is d*n*n = d*n^2
    cost_matmul1 = d * n**2
    
    # Cost of (...) @ X.T: (d,n) @ (n,d) -> (d,d). Cost is d*n*d = d^2*n
    cost_matmul2 = d**2 * n
    
    # The total cost is the sum of these major components.
    total_cost = cost_xtx + cost_inv + cost_gx + cost_matmul1 + cost_matmul2

    print("The minimum cost is achieved by avoiding the d^2 x d^2 inversion and instead ")
    print("inverting an n x n matrix. The primary costs are from matrix multiplications.")
    print("\nThe breakdown of the computational cost (in FLOPS) is as follows:")
    print(f"1. Cost of (X.T @ X)      : {d} * {n}^2 = {cost_xtx}")
    print(f"2. Cost of (n x n) inverse : {n}^3 = {cost_inv}")
    print(f"3. Cost of (G @ X)         : {d}^2 * {n} = {cost_gx}")
    print(f"4. Cost of 1st matmul      : {d} * {n}^2 = {cost_matmul1}")
    print(f"5. Cost of 2nd matmul      : {d}^2 * {n} = {cost_matmul2}")
    
    print("\nThe final equation for the total cost is:")
    print(f"Cost = ({d}*{n}^2) + ({n}^3) + ({d}^2*{n}) + ({d}*{n}^2) + ({d}^2*{n})")
    print(f"     = {cost_xtx} + {cost_inv} + {cost_gx} + {cost_matmul1} + {cost_matmul2}")
    print(f"     = {total_cost} FLOPS")

    print(f"\nSince n < d, the complexity is dominated by the d^2*n terms.")
    print("Therefore, the minimum achievable curvature cost is of the order O(d^2 * n).")


# Example values where d > n
d_example = 100
n_example = 20

calculate_curvature_cost(d_example, n_example)