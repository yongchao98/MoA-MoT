import itertools

def is_2_indistinguishable(s1, s2):
    """
    Checks if two strings are indistinguishable by any 2-state DFA.
    This is true iff they have the same last character and same parity of counts
    for each symbol.
    """
    if s1[-1] != s2[-1]:
        return False
    if s1.count('0') % 2 != s2.count('0') % 2:
        return False
    if s1.count('1') % 2 != s2.count('1') % 2:
        return False
    return True

def is_3_distinguishable(s1, s2):
    """
    Checks if two strings can be distinguished by at least one 3-state DFA.
    It iterates through all possible 3-state DFAs to find one that separates the strings.
    """
    num_states = 3
    # A DFA is defined by its transition function. For a 3-state, 2-symbol DFA,
    # there are 3*2=6 transitions to define. Each can go to one of 3 states.
    # So there are 3^6 = 729 possible DFAs.
    num_dfas = num_states ** (num_states * 2)

    for i in range(num_dfas):
        # Build the transition table for the i-th DFA
        temp = i
        transitions = {}
        for state in range(num_states):
            transitions[state] = {}
            for symbol_val in range(2):
                symbol = str(symbol_val)
                transitions[state][symbol] = temp % num_states
                temp //= num_states
        
        # Simulate the DFA on s1
        current_state_s1 = 0
        for char in s1:
            current_state_s1 = transitions[current_state_s1][char]

        # Simulate the DFA on s2
        current_state_s2 = 0
        for char in s2:
            current_state_s2 = transitions[current_state_s2][char]
        
        if current_state_s1 != current_state_s2:
            # Found a DFA that distinguishes them
            return True

    # No 3-state DFA could distinguish them
    return False

def find_minimum_n():
    """
    Finds the minimum length n where a pair of strings is 2-indistinguishable
    but 3-distinguishable.
    """
    n = 1
    while True:
        # Generate all binary strings of length n
        strings = ["".join(seq) for seq in itertools.product("01", repeat=n)]
        
        # Check all unique pairs of strings
        for s1, s2 in itertools.combinations(strings, 2):
            if is_2_indistinguishable(s1, s2):
                # Found a 2-indistinguishable pair. Now check if it's 3-distinguishable.
                if is_3_distinguishable(s1, s2):
                    # This is the minimum n that satisfies the conditions
                    return n
        n += 1

if __name__ == "__main__":
    min_n = find_minimum_n()
    print(min_n)