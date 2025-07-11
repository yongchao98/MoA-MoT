def solve():
    """
    Calculates the next three terms of the sequence starting from NZX.
    """
    # Initial state from the last given term: NZX
    # N -> 13, Z -> 25, X -> 23
    r_i, r_p, r_q = 13, 25, 23
    
    # Header for the output
    print("The last given term is NZX, which corresponds to (i=13, p=25, q=23).")
    print("Calculating the next three terms:")
    print("-" * 30)

    for step in range(1, 4):
        p_old = r_p
        
        # Calculate next state candidates
        q_new = (r_p + r_q) % 26
        p_new = (r_i + q_new) % 26
        
        # Letters for printing calculations
        char_i, char_p, char_q = chr(r_i+65), chr(r_p+65), chr(r_q+65)
        
        print(f"Step {step}: Starting from {char_i}{char_p}{char_q} (i={r_i}, p={r_p}, q={r_q})")
        print(f"  q_new = ({r_p} + {r_q}) % 26 = {q_new}")
        print(f"  p_new = ({r_i} + {q_new}) % 26 = {p_new}")

        if p_new < p_old:
            print(f"  Condition p_new({p_new}) < p_old({p_old}) is TRUE. Advancing to the next group.")
            r_i = r_i + 1
            # Recalculate p and q based on new i
            p_next = (r_i + q_new) % 26
            q_next = (p_next + q_new) % 26
            print(f"  New i = {r_i}")
            print(f"  New p = ({r_i} + {q_new}) % 26 = {p_next}")
            print(f"  New q = ({p_next} + {q_new}) % 26 = {q_next}")
            r_p = p_next
            r_q = q_next
        else:
            print(f"  Condition p_new({p_new}) < p_old({p_old}) is FALSE. Staying in the same group.")
            r_p = p_new
            r_q = q_new
            print(f"  New i = {r_i}")
            print(f"  New p = {r_p}")
            print(f"  New q = {r_q}")

        next_term_str = f"{chr(r_i+65)}{chr(r_p+65)}{chr(r_q+65)}"
        print(f"  Resulting term: {next_term_str}\n")
    
    print("The three capital letters that come next are: OKG, PFV, PPA")


solve()
<<<OKG, PFV, PPA>>>