def check_if_pair_can_form_infinite_group(n1, n2):
    """
    Checks if the pair (n1, n2) can form an infinite group GG_{n1, n2}(r) for some r > 0.
    Based on mathematical analysis, this is true if and only if both n1 and n2 are >= 2.
    """
    return n1 >= 2 and n2 >= 2

def find_min_N():
    """
    Finds the minimum positive integer N for which the set S(N) is non-empty.
    """
    N = 1
    while True:
        print(f"Checking N = {N}")
        is_S_N_non_empty = False
        found_pair = None
        # Iterate through all pairs (n1, n2) such that n1 <= N and n2 <= N
        for n1 in range(1, N + 1):
            for n2 in range(1, N + 1):
                # Check if this pair can form an infinite group
                if check_if_pair_can_form_infinite_group(n1, n2):
                    is_S_N_non_empty = True
                    found_pair = (n1, n2)
                    break
            if is_S_N_non_empty:
                break
        
        if is_S_N_non_empty:
            print(f"For N = {N}, S(N) is not empty.")
            print(f"A valid pair (n1, n2) is {found_pair}, since n1 = {found_pair[0]} >= 2 and n2 = {found_pair[1]} >= 2.")
            print(f"Also, n1 = {found_pair[0]} <= {N} and n2 = {found_pair[1]} <= {N}.")
            print(f"Therefore, the minimum N is {N}.")
            return N
        else:
            print(f"For N = {N}, S(N) is empty. No pair (n1, n2) with n1, n2 <= {N} satisfies the condition (n1 >= 2 and n2 >= 2).")
            N += 1

# Execute the search and print the final answer
min_N = find_min_N()