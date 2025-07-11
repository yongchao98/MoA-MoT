def count_reflection_types(k_values):
    """Counts reflections of type 1 (even k) and type 3/4 (odd k)."""
    count_even = 0
    count_odd = 0
    for k in k_values:
        if k % 2 == 0:
            count_even += 1
        else:
            count_odd += 1
    return count_even, count_odd

def solve():
    """
    Calculates the number of possible light ray paths based on reflection counts.
    """
    # Given number of reflections on each mirror
    n1_req = 2  # G1 reflections (vertical, even k)
    n3_req = 1  # G3 reflections (vertical, odd k)
    n2_req = 1  # G2 reflections (horizontal, even k)
    n4_req = 1  # G4 reflections (horizontal, odd k)

    # Total number of vertical and horizontal reflections
    p_mag = n1_req + n3_req
    q_mag = n2_req + n4_req

    possible_p = [p_mag, -p_mag]
    possible_q = [q_mag, -q_mag]

    valid_p_count = 0
    valid_q_count = 0

    print("Checking possible values for horizontal cell index p...")
    for p in possible_p:
        if p > 0:
            # Walls from k=1 to p
            k_vals_p = range(1, p + 1)
        else: # p < 0
            # Walls from k=0 down to p+1
            k_vals_p = range(0, p - 1, -1)
        
        n1_calc, n3_calc = count_reflection_types(k_vals_p)
        
        if n1_calc == n1_req and n3_calc == n3_req:
            valid_p_count += 1
            print(f"p = {p} is a valid choice (n1={n1_calc}, n3={n3_calc}).")
        else:
            print(f"p = {p} is not a valid choice (n1={n1_calc}, n3={n3_calc}).")

    print("\nChecking possible values for vertical cell index q...")
    for q in possible_q:
        if q > 0:
            # Walls from k=1 to q
            k_vals_q = range(1, q + 1)
        else: # q < 0
            # Walls from k=0 down to q+1
            k_vals_q = range(0, q - 1, -1)

        n2_calc, n4_calc = count_reflection_types(k_vals_q)
        
        if n2_calc == n2_req and n4_calc == n4_req:
            valid_q_count += 1
            print(f"q = {q} is a valid choice (n2={n2_calc}, n4={n4_calc}).")
        else:
            print(f"q = {q} is not a valid choice (n2={n2_calc}, n4={n4_calc}).")
    
    total_ways = valid_p_count * valid_q_count
    
    print("\n--- Calculation Summary ---")
    print(f"Number of valid p values: {valid_p_count}")
    print(f"Number of valid q values: {valid_q_count}")
    print(f"Total number of ways = {valid_p_count} * {valid_q_count} = {total_ways}")
    
solve()
<<<2>>>