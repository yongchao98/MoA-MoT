import sys

def solve():
    """
    Solves the remainder maximization problem using recursion with memoization.
    """
    # Increase recursion limit for larger n, if needed.
    # sys.setrecursionlimit(2000)

    # Example Input:
    # line 1: n (number of elements) and initial x
    # line 2: n numbers in the list a
    #
    # Example 1:
    # 2 10
    # 3 7
    # Expected result: permutation (3, 7) giving final x = 1
    # 10 % 3 = 1; 1 % 7 = 1;
    #
    # Example 2:
    # 3 25
    # 10 4 18
    # Expected result: permutation (18, 10, 4) giving final x = 3
    # 25 % 18 = 7; 7 % 10 = 7; 7 % 4 = 3;

    try:
        n_str, x_str = input().split()
        n = int(n_str)
        initial_x = int(x_str)
        a = list(map(int, input().split()))
    except (IOError, ValueError):
        # Default values if input fails or is empty
        n = 3
        initial_x = 25
        a = [10, 4, 18]
        print("Using default example: n=3, x=25, a=[10, 4, 18]")


    memo = {}

    def find_best_permutation(current_x, mask):
        """
        Recursively finds the permutation that maximizes the final remainder.

        :param current_x: The current value of x.
        :param mask: A bitmask representing the set of available numbers in 'a'.
        :return: A tuple (max_remainder, optimal_permutation).
        """
        if mask == 0:
            return current_x, []

        if (current_x, mask) in memo:
            return memo[(current_x, mask)]

        best_val = -1
        best_path = []

        for i in range(n):
            # Check if the i-th number is in the current subset
            if (mask >> i) & 1:
                # Recursively call for the subproblem without the i-th number
                val, path = find_best_permutation(current_x, mask ^ (1 << i))
                
                # Apply the current number a[i] to the result of the subproblem
                current_permutation_val = val % a[i]
                
                # We need to explore all paths from this state
                final_val, final_path_suffix = find_best_permutation(current_permutation_val, mask ^ (1 << i))
                
                if final_val > best_val:
                    best_val = final_val
                    best_path = [a[i]] + final_path_suffix


        memo[(current_x, mask)] = (best_val, best_path)
        return best_val, best_path
        
    # The problem statement can be interpreted in two ways for the DP state transition.
    # 1. State(x, mask): apply one number `a[i]`, then solve for state(x % a[i], mask \ {i}).
    # 2. State(mask): solve for subproblem `mask \ {i}`, get a value `v`, then compute `v % a[i]`.
    # Let's try the second, simpler DP formulation, which is standard for this class of problems.
    
    memo_simple = {}
    def find_best_permutation_simple(mask):
        """
        A simpler DP formulation. State depends only on the set of numbers used.
        The initial `x` is implicit.
        
        :param mask: Bitmask of numbers from `a` to use.
        :return: A tuple (max_remainder, optimal_permutation).
        """
        if mask == 0:
            return initial_x, []
            
        if mask in memo_simple:
            return memo_simple[mask]

        best_val = -1
        best_path = []
        
        for i in range(n):
            if (mask >> i) & 1:
                # Solve subproblem for the set without a[i]
                prev_val, prev_path = find_best_permutation_simple(mask ^ (1 << i))
                
                # Apply a[i] as the last operation for this permutation
                current_val = prev_val % a[i]
                
                # Update best result if this path is better
                if current_val > best_val:
                    best_val = current_val
                    best_path = prev_path + [a[i]]

        memo_simple[mask] = (best_val, best_path)
        return best_val, best_path

    # We use the simpler and more standard DP formulation.
    final_x, p = find_best_permutation_simple((1 << n) - 1)
    
    # Print the equation step-by-step
    equation = []
    current_val = initial_x
    equation.append(f"{current_val}")
    
    for val_a in p:
        prev_val = current_val
        current_val %= val_a
        equation.append(f"% {val_a} = {current_val}")

    print(" ".join(equation))

solve()