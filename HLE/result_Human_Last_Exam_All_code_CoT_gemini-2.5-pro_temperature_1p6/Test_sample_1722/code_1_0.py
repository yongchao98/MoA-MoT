import itertools

def count_substring(s, sub):
    """Counts non-overlapping occurrences of a substring."""
    return s.count(sub)

def are_2_indistinguishable(u, v):
    """
    Checks if two strings u and v are indistinguishable by any 2-state DFA.
    This is true iff they have the same start/end symbols and same parity of '01' and '10' counts.
    """
    if len(u) != len(v) or u == v:
        return False
    if u[0] != v[0] or u[-1] != v[-1]:
        return False
    if (count_substring(u, '01') % 2) != (count_substring(v, '01') % 2):
        return False
    if (count_substring(u, '10') % 2) != (count_substring(v, '10') % 2):
        return False
    return True

def are_3_distinguishable(u, v):
    """
    Checks if two strings u and v are distinguishable by some 3-state DFA.
    This is done by iterating through all possible 3-state DFAs.
    """
    num_states = 3
    num_symbols = 2 # Alphabet {0, 1}
    start_state = 0 # Conventionally, the first state is the start state.

    # A DFA's transition function can be represented as a list of lists.
    # We generate all possible transition functions.
    # Total DFAs = (num_states)^(num_states * num_symbols) = 3^(3*2) = 729.
    
    # This represents all possible transition rules for a single state: e.g., (state_on_0, state_on_1)
    state_transitions = list(itertools.product(range(num_states), repeat=num_symbols))
    
    # This generates all combinations of transition rules for all states.
    for dfa_config in itertools.product(state_transitions, repeat=num_states):
        # A specific 3-state DFA configuration
        
        # Simulate string u
        current_state_u = start_state
        for symbol in u:
            # dfa_config[current_state] is the tuple of transitions for that state
            # e.g., (next_on_0, next_on_1)
            current_state_u = dfa_config[current_state_u][int(symbol)]
            
        # Simulate string v
        current_state_v = start_state
        for symbol in v:
            current_state_v = dfa_config[current_state_v][int(symbol)]

        # If final states differ, they are distinguishable
        if current_state_u != current_state_v:
            return True
            
    # If no DFA could distinguish them
    return False

def find_minimum_length():
    """
    Finds the minimum length n for which there exists a pair of strings that are
    2-indistinguishable but 3-distinguishable.
    """
    n = 1
    while True:
        print(f"Searching with length n = {n}...")
        
        # Generate all binary strings of length n
        strings = ["".join(p) for p in itertools.product('01', repeat=n)]
        
        # Check all unique pairs of strings
        for u, v in itertools.combinations(strings, 2):
            if are_2_indistinguishable(u, v):
                print(f"  Found 2-indistinguishable pair: ({u}, {v})")
                if are_3_distinguishable(u, v):
                    print(f"  Pair ({u}, {v}) is also 3-distinguishable.")
                    print("\nThis is the first length that satisfies all conditions.")
                    print("Therefore, the minimum length of the hallway n is:")
                    # This fulfills the "output each number in the final equation" requirement
                    # by printing the final resulting number.
                    print(n)
                    return n
                else:
                    print(f"  But pair ({u}, {v}) is NOT 3-distinguishable.")

        n += 1

# Execute the search
find_minimum_length()
