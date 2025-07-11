def solve():
    """
    Calculates the sequence S4 based on the deduced rule and prints the equation for the 42nd term.
    """
    s = {1: 1, 2: 1}
    n_max = 42

    for n in range(3, n_max + 1):
        try:
            val_n_minus_1 = s[n - 1]
            val_n_minus_2 = s[n - 2]
            
            index1 = n - val_n_minus_1
            index2 = n - 1 - val_n_minus_2
            
            if index1 < 1 or index2 < 1:
                # This case might happen with some rules, handle gracefully.
                # For this rule, it shouldn't be an issue with n > 2.
                s[n] = -1 # Error indicator
                continue

            s[n] = s[index1] + s[index2]

        except KeyError:
            # This handles cases where an index is not yet computed.
            s[n] = -1 # Error indicator
            continue

    s41 = s.get(n_max - 1, 'N/A')
    s40 = s.get(n_max - 2, 'N/A')
    s42 = s.get(n_max, 'N/A')
    
    val_s41 = s41
    val_s40 = s40
    
    idx1 = n_max - val_s41
    idx2 = n_max - 1 - val_s40
    
    val1 = s.get(idx1, 'N/A')
    val2 = s.get(idx2, 'N/A')

    print(f"s[{n_max}] = s[{n_max} - s[{n_max - 1}]] + s[{n_max - 1} - s[{n_max - 2}]] = s[{n_max} - {s41}] + s[{n_max - 1} - {s40}] = s[{idx1}] + s[{idx2}] = {val1} + {val2} = {s42}")

solve()
<<<16>>>