def solve():
    """
    This function generates the sequence S4 based on a deduced recurrence relation.
    """
    s = {0: 0, 1: 1, 2: 1} # Using a dictionary to store sequence values, with a seed for s[0]
    n_terms = 42
    
    # Generate the sequence from the 3rd term up to n_terms
    for n in range(3, n_terms + 1):
        try:
            # The recurrence relation for S4
            # s[n] = s[n - s[n-1]] + s[s[n-2] - 1]
            term1_idx = n - s[n-1]
            term2_idx = s[n-2] - 1

            s[n] = s[term1_idx] + s[term2_idx]
        except KeyError as e:
            print(f"Error: index {e} not found in sequence while calculating s[{n}].")
            return
    
    # Print the resulting sequence
    result = [s[i] for i in range(1, n_terms + 1)]
    print("Generated Sequence S4:")
    print(", ".join(map(str, result)))
    
    # Print the final rule
    print("\nThe rule for S4 is:")
    print("s[1] = 1")
    print("s[2] = 1")
    print("R(s[n]) = s[n - s[n-1]] + s[s[n-2] - 1] (with s[0]=0)")

solve()
<<<R(s[n]) = s[n - s[n-1]] + s[s[n-2] - 1]>>>