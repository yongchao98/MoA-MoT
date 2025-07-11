import math

def solve_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polyhedral approximation B with edge lengths a.

    Args:
        n (int): The number of sides of the polygon B.
        a (list of float): The edge lengths [a_1, a_2, ..., a_n] of B.
    """
    if len(a) != n:
        print("Error: The number of edge lengths must be equal to n.")
        return

    max_x = 0

    if n % 2 == 1:
        # n is odd, solution is unique
        x = [0] * n
        for i in range(n):
            current_sum = 0
            for j in range(n):
                idx = (i + j) % n
                current_sum += ((-1)**j) * a[idx]
            x[i] = 0.5 * current_sum
        
        # Check for validity (all x_i must be non-negative)
        if all(val >= -1e-9 for val in x):
            max_x = max(x)
        else:
            print("The given edge lengths do not form a valid circumscribed polygon.")
            return
            
    else:
        # n is even, solution is not unique.
        # Check if a solution can exist.
        sum_odd_edges = sum(a[1::2])
        sum_even_edges = sum(a[0::2])
        if not math.isclose(sum_odd_edges, sum_even_edges, rel_tol=1e-9):
            print("Error: For even n, the sum of odd and even edge lengths must be equal.")
            print(f"Sum of even-indexed edges: {sum_even_edges}")
            print(f"Sum of odd-indexed edges: {sum_odd_edges}")
            return
            
        global_max_x = 0
        # Find the max value by checking boundary solutions (where one x_k = 0).
        for k in range(n):
            x = [0] * n
            x[k] = 0
            
            # Propagate the solution from x[k]
            # using the relation x_{i+1} = a_i - x_i
            for i in range(n - 1):
                current_idx = (k + i) % n
                next_idx = (k + i + 1) % n
                x[next_idx] = a[current_idx] - x[current_idx]

            # Check if this boundary solution is valid (all x_i >= 0)
            if all(val >= -1e-9 for val in x):
                current_max_x = max(x)
                if current_max_x > global_max_x:
                    global_max_x = current_max_x
        
        max_x = global_max_x

    sin_pi_over_n = math.sin(math.pi / n)
    result = max_x * sin_pi_over_n
    
    print("Final equation: max_dist = max(x_i) * sin(pi/n)")
    print(f"max(x_i) = {max_x}")
    print(f"sin(pi/n) = sin({math.pi/n:.4f}) = {sin_pi_over_n}")
    print(f"Largest possible Hausdorff distance = {max_x} * {sin_pi_over_n} = {result}")
    
    # Return the final answer in the requested format
    print(f"\n<<<{result}>>>")

# --- Example Usage ---
# We will solve for a case with n=4, which is even and more complex.
# The condition is a_0 + a_2 = a_1 + a_3.
# Let a = [10, 8, 12, 14]. Then 10 + 12 = 22 and 8 + 14 = 22. Condition holds.
n_example = 4
a_example = [10.0, 8.0, 12.0, 14.0]

print(f"Solving for n={n_example} and edge lengths a={a_example}\n")
solve_hausdorff_distance(n_example, a_example)
