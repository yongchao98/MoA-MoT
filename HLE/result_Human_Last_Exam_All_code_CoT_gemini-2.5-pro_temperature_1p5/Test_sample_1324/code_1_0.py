import math

def simulate_step(chocolates):
    """Simulates one minute of the chocolate passing game."""
    n = len(chocolates)
    new_chocolates = [0] * n
    for i in range(n):
        left_neighbor_idx = (i - 1 + n) % n
        # Calculate the sum from the person and their left neighbor
        s = chocolates[i] + chocolates[left_neighbor_idx]
        
        # Each person passes half, so they receive half from the left and keep half of their own
        intermediate_value = s // 2
        
        # If the result is odd, they take an extra chocolate
        if intermediate_value % 2 != 0:
            new_chocolates[i] = intermediate_value + 1
        else:
            new_chocolates[i] = intermediate_value
            
    return new_chocolates

def analyze_simulation(n, initial_chocolates, steps):
    """Runs a simulation and analyzes the statements."""
    chocolates = list(initial_chocolates)
    history = []

    print(f"Initial state C(0): {chocolates}")
    
    l_i, h_i = min(chocolates), max(chocolates)
    d_i = h_i - l_i
    history.append({'l': l_i, 'h': h_i, 'd': d_i})
    print(f"i=0: l={l_i}, h={h_i}, d={d_i}, C={chocolates}")
    
    # Run simulation
    for i in range(1, steps):
        chocolates = simulate_step(chocolates)
        l_i, h_i = min(chocolates), max(chocolates)
        d_i = h_i - l_i
        history.append({'l': l_i, 'h': h_i, 'd': d_i})
        print(f"i={i}: l={l_i}, h={h_i}, d={d_i}, C={chocolates}")

    # --- Statement Analysis ---
    print("\n--- Analyzing Statements ---")
    
    # Statement 1 check for i=0
    s1_holds_i0 = False
    for m in range(1, n):
        if m < len(history):
            if history[m]['d'] < history[0]['d']:
                s1_holds_i0 = True
                print(f"Statement 1 check (i=0): TRUE. For m={m}, d({m})={history[m]['d']} < d(0)={history[0]['d']}")
                break
    if not s1_holds_i0:
        print(f"Statement 1 check (i=0): FALSE. No m in [1, {n-1}] found where d(m) < d(0).")

    # Statement 2 check
    s2_holds = False
    for i in range(1, len(history)):
        if history[i]['l'] < history[i-1]['l']:
            s2_holds = True
            print(f"Statement 2 check: TRUE. At i={i}, l({i})={history[i]['l']} < l({i-1})={history[i-1]['l']}.")
            break
    if not s2_holds:
        print("Statement 2 check: FALSE. The minimum number of chocolates never decreased.")

    # Statement 3 check for i=0
    s3_holds_i0 = False
    for m in range(1, n):
        if m < len(history):
            if history[m]['l'] > history[0]['l']:
                s3_holds_i0 = True
                print(f"Statement 3 check (i=0): TRUE. For m={m}, l({m})={history[m]['l']} > l(0)={history[0]['l']}")
                break
    if not s3_holds_i0:
        print(f"Statement 3 check (i=0): FALSE. No m in [1, {n-1}] found where l(m) > l(0).")


if __name__ == '__main__':
    print("--- Case 1: Non-stable initial state ---")
    n1 = 4
    initial1 = [10, 20, 6, 14]
    analyze_simulation(n1, initial1, steps=5)

    print("\n" + "="*50 + "\n")

    print("--- Case 2: Stable initial state ---")
    n2 = 4
    initial2 = [10, 10, 10, 10]
    analyze_simulation(n2, initial2, steps=5)

    print("\n" + "="*50 + "\n")
    print("Final Conclusion:")
    print("Statement 2 is provably false as the minimum l(i) is non-decreasing.")
    print("Statements 1 and 3 are phrased with 'For any i >= 0', which implies they must hold for all possible configurations, including stable ones.")
    print("As shown in Case 2 (stable state), both Statements 1 and 3 fail.")
    print("Therefore, none of the statements are true.")
