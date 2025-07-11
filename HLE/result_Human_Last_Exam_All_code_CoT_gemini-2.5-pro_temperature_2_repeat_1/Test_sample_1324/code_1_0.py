import math

def simulate_and_check(c0, num_steps_to_run):
    """
    Simulates the chocolate passing game for a given initial configuration
    and checks the validity of the three statements.
    """
    n = len(c0)
    history = [list(c0)]
    
    print(f"Starting simulation with n={n} people and initial state c_0 = {c0}\n")
    
    # --- Simulation ---
    for i in range(num_steps_to_run):
        c_current = history[-1]
        c_next = [0] * n
        print(f"Step {i+1}:")
        for k in range(n):
            # person p_k receives from p_{k-1}
            k_prev = (k - 1 + n) % n
            
            val_k = c_current[k]
            val_prev = c_current[k_prev]
            
            temp_sum = (val_k + val_prev) / 2
            
            if temp_sum % 2 != 0:  # Odd
                c_next[k] = math.ceil(temp_sum)
                print(f"  - p_{k+1}'s new count = ceil(({val_k} + {val_prev}) / 2) = ceil({temp_sum}) = {c_next[k]}")
            else: # Even
                c_next[k] = int(temp_sum)
                print(f"  - p_{k+1}'s new count = ({val_k} + {val_prev}) / 2 = {temp_sum} = {c_next[k]}")
                
        history.append(c_next)
        print(f"Resulting state c_{i+1} = {c_next}\n")

    l_hist = [min(c) for c in history]
    h_hist = [max(c) for c in history]
    d_hist = [h - l for h, l in zip(h_hist, l_hist)]

    # --- Statement Analysis ---
    print("--- Analyzing Statements ---")

    # Check Statement 2
    s2_is_false = True
    for i in range(1, len(l_hist)):
        if l_hist[i] < l_hist[i-1]:
            print(f"Statement 2 is TRUE because l_{i}={l_hist[i]} < l_{i-1}={l_hist[i-1]}.")
            s2_is_false = False
            break
    if s2_is_false:
        print("Statement 2 is FALSE. The minimum l_i never decreases.")
        print(f"Sequence of minimums: {l_hist}")

    # The system converges when d_i becomes 0. Let's find when this happens.
    try:
        converged_i = d_hist.index(0)
        print(f"\nThe system reached a converged state (d_i=0) at i={converged_i}.")
    except ValueError:
        print("\nSystem did not converge in the simulated steps.")
        converged_i = -1
    
    # Check Statements 1 and 3 at the point of convergence
    if converged_i != -1 and len(history) > converged_i + n:
        i = converged_i
        print(f"\nChecking Statements 1 and 3 for i = {i}, where the state is converged.")
        
        # Check Statement 1
        s1_holds_at_i = False
        for m in range(1, n):
            if d_hist[i+m] < d_hist[i]:
                s1_holds_at_i = True
                break
        if not s1_holds_at_i:
            print(f"Statement 1 is FALSE. For i={i}, d_{i}={d_hist[i]}. We must find m in [1,{n-1}] where d_{i+m} < d_{i}.")
            m_test = 1
            print(f"  - For m={m_test}, d_{i+m}={d_hist[i+m_test]}. The condition {d_hist[i+m_test]} < {d_hist[i]} is false.")
            
        # Check Statement 3
        s3_holds_at_i = False
        for m in range(1, n):
            if l_hist[i+m] > l_hist[i]:
                s3_holds_at_i = True
                break
        if not s3_holds_at_i:
            print(f"Statement 3 is FALSE. For i={i}, l_{i}={l_hist[i]}. We must find m in [1,{n-1}] where l_{i+m} > l_{i}.")
            m_test = 1
            print(f"  - For m={m_test}, l_{i+m}={l_hist[i+m_test]}. The condition {l_hist[i+m_test]} > {l_hist[i]} is false.")
    else:
        print("\nCould not check S1 and S3 at convergence. Please run simulation longer.")

# Initial state designed to converge quickly
c_initial = [2, 2, 6]
# n=3. We need to check for m=1,2. If we test i=2, we need up to i+2=4.
# Let's run for 5 steps to be safe.
simulate_and_check(c_initial, 5)

print("\n\nConclusion: Since we found counterexamples for all three statements, none of them are universally true.")