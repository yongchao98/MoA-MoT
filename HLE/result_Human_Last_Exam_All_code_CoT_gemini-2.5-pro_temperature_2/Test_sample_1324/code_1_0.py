import math

def simulate_chocolates(initial_chocolates, num_steps):
    """
    Simulates the chocolate passing game and analyzes the statements.
    """
    n = len(initial_chocolates)
    chocolates = list(initial_chocolates)
    
    history_l = []
    history_d = []

    print(f"Starting simulation with n={n} and initial state: {chocolates}\n")

    for i in range(num_steps + 1):
        l = min(chocolates)
        h = max(chocolates)
        d = h - l
        history_l.append(l)
        history_d.append(d)
        
        print(f"Minute i={i}: Chocolates={chocolates}, l^{i}={l}, d^{i}={d}")

        # Calculate the next state
        next_chocolates = [0] * n
        for k in range(n):
            # Person k receives from person k-1 (with wrap-around)
            left_neighbor_idx = (k - 1 + n) % n
            
            # Sum of chocolates from self (half) and neighbor (half)
            temp_sum = chocolates[k] + chocolates[left_neighbor_idx]
            new_val = temp_sum / 2
            
            # Add one if the result is odd
            if new_val % 2 != 0:
                new_val += 1
            
            next_chocolates[k] = int(new_val)
        chocolates = next_chocolates

    print("\n--- Analyzing Statements ---")
    
    # Analysis of Statement 1
    # Check if for any i, d[i+m] < d[i] for some m < n
    s1_holds = True
    for i in range(num_steps - n + 1):
        found_m = False
        for m in range(1, n):
            if (i + m) < len(history_d) and history_d[i+m] < history_d[i]:
                found_m = True
                break
        if not found_m:
            s1_holds = False
            print(f"Statement 1 is FALSE. Counterexample found at i={i}.")
            print(f"For i={i}, d^{i}={history_d[i]}. We need d^{i+m} < {history_d[i]} for some m in [1, {n-1}].")
            for m in range(1, n):
                 if (i + m) < len(history_d):
                     print(f"  m={m}: d^{i+m}={history_d[i+m]}. The condition {history_d[i+m]} < {history_d[i]} is false.")
            break
    if s1_holds:
        print("Statement 1 appears to hold for the simulated steps (but may be false in general).")

    # Analysis of Statement 2
    # Check if l[i] < l[i-1] for any i > 0
    s2_holds = False
    for i in range(1, len(history_l)):
        if history_l[i] < history_l[i-1]:
            s2_holds = True
            print(f"Statement 2 is TRUE. Found case at i={i}: l^{i}={history_l[i]} < l^{i-1}={history_l[i-1]}")
            break
    if not s2_holds:
        print(f"Statement 2 is FALSE. The minimum 'l' is non-decreasing in the simulation: {history_l}")

    # Analysis of Statement 3
    # Check if for any i, l[i+m] > l[i] for some m < n
    s3_holds = True
    for i in range(num_steps - n + 1):
        found_m = False
        for m in range(1, n):
            if (i + m) < len(history_l) and history_l[i+m] > history_l[i]:
                found_m = True
                break
        if not found_m:
            s3_holds = False
            print(f"Statement 3 is FALSE. Counterexample found at i={i}.")
            print(f"For i={i}, l^{i}={history_l[i]}. We need l^{i+m} > {history_l[i]} for some m in [1, {n-1}].")
            for m in range(1, n):
                 if (i + m) < len(history_l):
                     print(f"  m={m}: l^{i+m}={history_l[i+m]}. The condition {history_l[i+m]} > {history_l[i]} is false.")
            break
    if s3_holds:
         print("Statement 3 appears to hold for the simulated steps (but may be false in general).")

if __name__ == '__main__':
    # Use an equilibrium state as a clear counterexample
    n_people = 4
    # All must start with an even number of chocolates. Let's use 10.
    initial = [10] * n_people 
    # Number of simulation steps. n_people steps are enough to check the conditions.
    steps = n_people 
    simulate_chocolates(initial, steps)
    print("\nConclusion: Based on the analysis, none of the statements are true.")
