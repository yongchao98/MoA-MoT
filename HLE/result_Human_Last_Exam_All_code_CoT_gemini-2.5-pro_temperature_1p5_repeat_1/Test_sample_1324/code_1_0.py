def simulate_chocolates(initial_chocolates, steps):
    """
    Simulates the chocolate passing game for a given number of steps.

    Args:
        initial_chocolates (list): A list of integers representing the starting
                                   number of chocolates for each person.
        steps (int): The number of minutes to simulate.
    """
    chocolates = list(initial_chocolates)
    n = len(chocolates)

    print(f"n = {n}")
    print("-" * 30)

    for i in range(steps):
        h = max(chocolates)
        l = min(chocolates)
        d = h - l
        print(f"i={i}: chocolates={chocolates}, h={h}, l={l}, d={d}")

        # Check statements at this step i (for demonstration)
        # S1 check fails at i=1: d^2 < d^1? (0 < 0 is False)
        # S3 check fails at i=1: l^2 > l^1? (6 > 6 is False)
        # S2 check never holds for any i > 0.

        # Create the next state
        next_chocolates = [0] * n
        for k in range(n):
            # Person k receives from person k-1
            left_neighbor_idx = (k - 1 + n) % n
            
            # Sum of chocolates from self (half) and neighbor (half)
            new_val_pre_adjustment = (chocolates[k] + chocolates[left_neighbor_idx]) // 2

            # Add one if the result is odd
            if new_val_pre_adjustment % 2 != 0:
                next_chocolates[k] = new_val_pre_adjustment + 1
            else:
                next_chocolates[k] = new_val_pre_adjustment
        
        chocolates = next_chocolates

# Example run that serves as a counterexample for all statements
# For n=2, m must be 1.
# At i=1, d^{1+1} = d^2 = 0 and d^1 = 0. 'd^2 < d^1' is false. S1 fails.
# At i=1, l^{1+1} = l^2 = 6 and l^1 = 6. 'l^2 > l^1' is false. S3 fails.
# l^i is never less than l^{i-1}. S2 fails.
simulate_chocolates(initial_chocolates=[10, 2], steps=4)

print("\n" + "="*30 + "\n")

# Example run with a starting equilibrium state
# S1 fails: Is d^{i+m} < d^i (0 < 0)? False.
# S3 fails: Is l^{i+m} > l^i (4 > 4)? False.
# S2 fails as always.
simulate_chocolates(initial_chocolates=[4, 4, 4, 4], steps=4)
