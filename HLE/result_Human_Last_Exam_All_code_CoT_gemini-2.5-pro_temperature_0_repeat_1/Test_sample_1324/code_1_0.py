import sys

def run_simulation_and_check_statements(initial_chocolates, n_steps_to_run):
    """
    Runs the chocolate passing simulation and checks the three statements.

    Args:
        initial_chocolates (list): The initial number of chocolates for each person.
        n_steps_to_run (int): The number of minutes to simulate.
    """
    n = len(initial_chocolates)
    if n <= 1:
        print("Number of people n must be greater than 1.")
        return

    # Store the history of chocolates, highs, lows, and differences
    history_c = [initial_chocolates]
    history_h = [max(initial_chocolates)]
    history_l = [min(initial_chocolates)]
    history_d = [history_h[0] - history_l[0]]

    print(f"Initial state c[0]: {history_c[0]}, l[0]: {history_l[0]}, d[0]: {history_d[0]}\n")

    # Run the simulation
    current_c = list(initial_chocolates)
    for i in range(n_steps_to_run):
        next_c = [0] * n
        for p_idx in range(n):
            # p_idx is the person receiving chocolates (person b in the formula)
            # prev_p_idx is the person passing chocolates (person a in the formula)
            prev_p_idx = (p_idx - 1 + n) % n
            
            # Calculate the sum
            temp_sum = (current_c[p_idx] + current_c[prev_p_idx]) / 2.0
            
            # If the sum is odd, add 1 to make it even
            if temp_sum % 2 != 0:
                next_c[p_idx] = int(temp_sum + 1)
            else:
                next_c[p_idx] = int(temp_sum)
        
        current_c = next_c
        history_c.append(current_c)
        history_h.append(max(current_c))
        history_l.append(min(current_c))
        history_d.append(history_h[-1] - history_l[-1])

    # --- Statement Analysis ---
    print("--- Analyzing Statements ---")
    
    # Statement 2: After the i-th minute for some i > 0, l[i] < l[i-1].
    s2_is_true = False
    print("\nAnalysis of Statement 2: l[i] < l[i-1] for some i > 0?")
    for i in range(1, len(history_l)):
        print(f"Checking i={i}: Is l[{i}] < l[{i-1}]? Is {history_l[i]} < {history_l[i-1]}?")
        if history_l[i] < history_l[i-1]:
            s2_is_true = True
            print("  -> True. Statement 2 holds for this i.")
        else:
            print("  -> False.")
    if not s2_is_true:
        print("Conclusion: Statement 2 was not found to be true in this simulation. (Theoretical proof shows it's always false).")

    # Statements 1 and 3 must hold for ANY i >= 0. We test for i=0.
    # If they fail for i=0 for our chosen initial state, the statements are false.
    
    # Statement 1: For any i >= 0, d[i+m] < d[i] where m < n.
    s1_holds_for_i0 = False
    print("\nAnalysis of Statement 1 at i=0: Is d[m] < d[0] for some m in [1, n-1]?")
    d0 = history_d[0]
    print(f"d[0] = {d0}")
    for m in range(1, n):
        if m < len(history_d):
            dm = history_d[m]
            print(f"Checking m={m}: Is d[{m}] < d[0]? Is {dm} < {d0}?")
            if dm < d0:
                s1_holds_for_i0 = True
                print(f"  -> True. Statement 1 holds for i=0 with m={m}.")
                break
            else:
                print("  -> False.")
    if not s1_holds_for_i0:
        print(f"Conclusion: For i=0, we could not find an m < {n} where d[m] < d[0].")
        print("Since Statement 1 must hold for *any* i, this counterexample proves it is FALSE.")

    # Statement 3: For any i >= 0, l[i+m] > l[i] where m < n.
    s3_holds_for_i0 = False
    print("\nAnalysis of Statement 3 at i=0: Is l[m] > l[0] for some m in [1, n-1]?")
    l0 = history_l[0]
    print(f"l[0] = {l0}")
    for m in range(1, n):
        if m < len(history_l):
            lm = history_l[m]
            print(f"Checking m={m}: Is l[{m}] > l[0]? Is {lm} > {l0}?")
            if lm > l0:
                s3_holds_for_i0 = True
                print(f"  -> True. Statement 3 holds for i=0 with m={m}.")
                break
            else:
                print("  -> False.")
    if not s3_holds_for_i0:
        print(f"Conclusion: For i=0, we could not find an m < {n} where l[m] > l[0].")
        print("Since Statement 3 must hold for *any* i, this counterexample proves it is FALSE.")

    print("\n--- Final Conclusion ---")
    print("Statement 2 is false because the minimum number of chocolates l[i] can never decrease.")
    print("Statements 1 and 3 are false because they must hold for any state, but they fail for an equilibrium state (as shown in this simulation), which is a valid configuration.")
    print("Therefore, none of the statements are true.")


if __name__ == '__main__':
    # We choose an initial state that is already in equilibrium.
    # This serves as a direct counterexample for Statements 1 and 3.
    # The number of people, n
    n_people = 5
    # The initial number of chocolates for each person (must be even)
    # This is an equilibrium state.
    initial_c = [20] * n_people
    
    # We need to run for at least n-1 steps to check the conditions for m < n.
    run_simulation_and_check_statements(initial_c, n_steps_to_run=n_people)
    
    # The final answer based on the analysis.
    print("\n<<<A>>>")
