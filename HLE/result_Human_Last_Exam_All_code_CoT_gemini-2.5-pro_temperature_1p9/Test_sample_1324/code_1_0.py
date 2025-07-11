import math

def simulate_chocolate_passing(initial_chocolates, total_steps):
    """
    Simulates the chocolate passing game for a given number of steps
    and demonstrates why the statements are false.
    """
    chocolates = list(initial_chocolates)
    n = len(chocolates)

    print(f"Starting with n={n} people.")
    print("-" * 40)

    for i in range(total_steps):
        h_i = max(chocolates)
        l_i = min(chocolates)
        d_i = h_i - l_i
        
        print(f"Step {i}:")
        print(f"  Chocolates: {chocolates}")
        print(f"  h^{i}={h_i}, l^{i}={l_i}, d^{i}={d_i}")

        if h_i == l_i:
            print("\nSystem has reached equilibrium.")
            print("Let's check the statements for this step (i={i}):")
            
            # Check Statement 1
            # We need d^{i+m} < d^i for some m < n. Here d^i=0.
            # d^{i+m} will also be 0. So 0 < 0 is false.
            print("\n- Statement 1 ('d^{i+m} < d^i') fails:")
            print(f"  At step {i}, d^{i}=0. For any future step, d will remain 0.")
            print(f"  The condition becomes 0 < 0, which is false.")

            # Check Statement 2
            # l^i should never be less than l^{i-1}
            print("\n- Statement 2 ('l^i < l^{i-1}') is false:")
            print("  As proven in the analysis, the minimum count l^i is non-decreasing.")
            
            # Check Statement 3
            # We need l^{i+m} > l^i for some m < n. Here l^i = C.
            # l^{i+m} will also be C. So C > C is false.
            print("\n- Statement 3 ('l^{i+m} > l^i') fails:")
            print(f"  At step {i}, l^{i}={l_i}. For any future step, l will remain {l_i}.")
            print(f"  The condition becomes {l_i} > {l_i}, which is false.")
            
            print("-" * 40)
            break

        print("\n  Calculating new chocolate distribution:")
        next_chocolates = [0] * n
        for k in range(n):
            # p_k gives to p_{k+1}, so p_{k+1} receives from p_k.
            # Using 0-based indexing for people p_0 to p_{n-1}
            # Person `k` gets chocolates from person `k-1`
            person_k_chocs = chocolates[k]
            neighbor_chocs = chocolates[(k - 1 + n) % n]
            
            temp_total = (person_k_chocs + neighbor_chocs) / 2
            
            if temp_total % 2 != 0 and math.floor(temp_total) == temp_total: 
                # temp_total is an odd integer
                new_chocs = int(temp_total + 1)
                print(f"    p_{k}: ({person_k_chocs} + {neighbor_chocs}) / 2 = {temp_total} (odd) -> gets {new_chocs} chocolates")
            else:
                # temp_total is even or not an integer (which is impossible here)
                new_chocs = int(temp_total)
                print(f"    p_{k}: ({person_k_chocs} + {neighbor_chocs}) / 2 = {temp_total} (even) -> gets {new_chocs} chocolates")
            
            next_chocolates[k] = new_chocs
        
        chocolates = next_chocolates
        print("-" * 40)

# Initial state that converges to equilibrium relatively quickly.
# n = 4, p_0 to p_3
initial_state = [10, 10, 6, 6]
simulate_chocolate_passing(initial_state, 10)