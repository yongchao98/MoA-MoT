def simulate_chocolates(initial_chocolates, steps):
    """
    Simulates the chocolate passing game for a given number of steps.
    """
    chocolates = list(initial_chocolates)
    n = len(chocolates)
    print(f"Let n = {n}")
    print("-" * 30)

    for i in range(steps):
        h_i = max(chocolates)
        l_i = min(chocolates)
        d_i = h_i - l_i
        
        print(f"Minute i = {i}")
        print(f"  Chocolates c^{i}: {chocolates}")
        print(f"  Max h^{i} = {h_i}, Min l^{i} = {l_i}, Diff d^{i} = {d_i}")

        if i > 0:
            # Check for Statement 2 counterexample (l_i >= l_{i-1})
            if l_i < l_prev:
                 print(f"  Statement 2 check: l^{i} < l^{i-1} ({l_i} < {l_prev}) - This should not happen.")
            else:
                 print(f"  Statement 2 check: l^{i} >= l^{i-1} ({l_i} >= {l_prev}) - True, so S2 is false.")


        if i == 1:
            # Check for Statement 1 counterexample with i=0, m=1
            print(f"Statement 1 check for i=0, m=1 (since m=1 < n={n}): Is d^1 < d^0?")
            print(f"  Is {d_i} < {d_0}? {d_i < d_0}")
            if not d_i < d_0:
                print("  This provides a counterexample for Statement 1.")

        print("-" * 30)
        
        # Save values for comparison in the next iteration
        l_prev = l_i
        if i == 0:
            d_0 = d_i

        next_chocolates = [0] * n
        for k in range(n):
            # Each person k receives from person k-1
            # Their new total is based on their current chocolates and their left neighbor's
            left_neighbor_idx = (k - 1 + n) % n
            
            # The sum is half of their own plus half of their neighbor's
            total_before_jar = (chocolates[k] + chocolates[left_neighbor_idx]) // 2

            # If the number is odd, take 1 from the jar to make it even
            if total_before_jar % 2 != 0:
                next_chocolates[k] = total_before_jar + 1
            else:
                next_chocolates[k] = total_before_jar
        
        chocolates = next_chocolates

        # If we reach a stable state, stop early and check Statement 3
        if h_i == l_i:
            print(f"System reached a stable state at i = {i}.")
            print(f"Checking Statement 3 for this state with l^{i}={l_i}.")
            # The state will now be constant
            l_i_plus_m = l_i 
            print("For any m > 0, l^(i+m) will be " + str(l_i_plus_m) + ".")
            print(f"Is there an m < n such that l^(i+m) > l^i (i.e., {l_i_plus_m} > {l_i})? This is false.")
            print("This stable state provides a counterexample for Statement 3.")
            break

# Main execution
# Case 1: Counterexample for Statement 1
print("Running simulation for Statement 1 counterexample...")
simulate_chocolates([10, 10, 14, 14], 3)

print("\n\n")

# Case 2: Counterexample for Statement 3
print("Running simulation for Statement 3 counterexample...")
simulate_chocolates([12, 12, 12], 2)