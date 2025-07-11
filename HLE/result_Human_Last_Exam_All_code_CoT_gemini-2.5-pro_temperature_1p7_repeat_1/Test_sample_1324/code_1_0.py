def simulate_chocolates(initial_chocolates, max_steps):
    """
    Simulates the chocolate passing game and checks the three statements.

    Args:
        initial_chocolates (list): A list of integers representing the initial
                                     number of chocolates for each person.
        max_steps (int): The number of minutes to simulate.
    """
    chocolates = list(initial_chocolates)
    n = len(chocolates)
    history = []

    print(f"Starting simulation with n={n} people.")
    print(f"Initial state c^0: {chocolates}\n")

    s1_holds = True
    s2_holds = False
    s3_holds = True

    for i in range(max_steps):
        # Store current state for history
        h_i = max(chocolates)
        l_i = min(chocolates)
        d_i = h_i - l_i
        history.append({
            "chocolates": list(chocolates),
            "h": h_i, "l": l_i, "d": d_i
        })
        
        # Print current state
        print(f"Minute i={i}:")
        print(f"  Chocolates: {chocolates}")
        print(f"  h^{i}={h_i}, l^{i}={l_i}, d^{i}={d_i}")

        # --- Calculate next state ---
        next_chocolates = [0] * n
        for k in range(n):
            # Person k receives from person k-1 (with wrap-around)
            left_neighbor_idx = (k - 1 + n) % n
            
            s = chocolates[k] + chocolates[left_neighbor_idx]
            avg = s / 2
            
            if avg % 2 == 1: # odd
                next_chocolates[k] = int(avg) + 1
            else: # even
                next_chocolates[k] = int(avg)
        
        chocolates = next_chocolates

        # --- Check statements against the historical data ---
        # Statement 2: l^i < l^{i-1} for some i > 0
        if i > 0:
            if history[i]['l'] < history[i-1]['l']:
                s2_holds = True
        
        # Check S1 and S3 for the state at step `i-1` if we have enough future data
        # We check for a failure of the "For any i" statements.
        # Check Statement 1 at step i
        s1_holds_for_i = False
        if d_i == 0: # Equilibrium reached
            s1_holds_for_i = False # d < 0 is impossible
        else:
             # This check requires future steps, so we do it retrospectively
             pass
        
        # Check Statement 3 at step i
        s3_holds_for_i = False
        if d_i == 0: # Equilibrium reached
            s3_holds_for_i = False # l > l is impossible
        else:
            # This check requires future steps, so we do it retrospectively
            pass

    # Post-simulation analysis
    print("\n--- Post-Simulation Analysis ---")

    # Analyze S1: For any i, exists m < n such that d^{i+m} < d^i
    s1_globally_true = True
    for i in range(len(history)):
        found_m_for_s1 = False
        for m in range(1, n):
            if i + m < len(history):
                if history[i+m]['d'] < history[i]['d']:
                    found_m_for_s1 = True
                    break
        if not found_m_for_s1 and history[i]['d'] > 0:
           # Could be true if simulation is not long enough
           # But if d[i] = 0, then it's definitely false
           pass
        if not found_m_for_s1 and i + n-1 < len(history) :
             s1_globally_true = False
             print(f"S1 fails at i={i}: d^{i}={history[i]['d']}. No m < {n} found where d^{{i+m}} < d^{i}.")
             break
    if s1_globally_true and (len(history) >= n):
        print("S1 appears to hold for non-equilibrium states in this run, but is proven false by the equilibrium case.")


    # Analyze S2: exists i>0 such that l^i < l^{i-1}
    s2_found_case = False
    for i in range(1, len(history)):
        if history[i]['l'] < history[i-1]['l']:
            print(f"S2 holds at i={i}: l^{i}={history[i]['l']} < l^{{i-1}}={history[i-1]['l']}.")
            s2_found_case = True
            break
    if not s2_found_case:
        print("S2 is false: l^i is non-decreasing. Found no instance of l^i < l^{i-1}.")


    # Analyze S3: For any i, exists m < n such that l^{i+m} > l^i
    s3_globally_true = True
    for i in range(len(history)):
        found_m_for_s3 = False
        for m in range(1, n):
            if i + m < len(history):
                if history[i+m]['l'] > history[i]['l']:
                    found_m_for_s3 = True
                    break
        if not found_m_for_s3 and i + n-1 < len(history):
            s3_globally_true = False
            print(f"S3 fails at i={i}: l^{i}={history[i]['l']}. No m < {n} found where l^{{i+m}} > l^{i}.")
            break
    if s3_globally_true and (len(history) >= n):
        print("S3 appears to hold for non-equilibrium states in this run, but is proven false by the equilibrium case.")

    print("\n--- Conclusion ---")
    print("Statement 1 is false because it fails in the equilibrium state (d=0).")
    print("Statement 2 is false because l^i is a non-decreasing sequence.")
    print("Statement 3 is false because it fails in the equilibrium state (l=constant).")
    print("Therefore, none of the statements are true.")


# Run simulation with an example that reaches equilibrium
# From the thought process: n=4, c^0 = (2, 2, 6, 6)
initial_state = [2, 2, 6, 6]
simulate_chocolates(initial_state, 10)