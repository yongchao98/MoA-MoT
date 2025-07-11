def solve_cellular_automaton():
    """
    Finds the intermediate step in a 2-generation cellular automaton evolution
    by testing all 256 possible elementary rules.
    """
    initial_state = '01101001'
    target_final_state = '10000111'

    for rule_num in range(256):
        # Calculate the intermediate state (generation t+1)
        intermediate_state = get_next_generation(initial_state, rule_num)
        
        # Calculate the next state from the intermediate one (generation t+2)
        calculated_final_state = get_next_generation(intermediate_state, rule_num)
        
        # Check if this rule produces the correct final state
        if calculated_final_state == target_final_state:
            # The problem states there's only one valid solution, so we print and exit.
            print(intermediate_state)
            return

def get_next_generation(state, rule_num):
    """
    Computes the next generation of a 1D cellular automaton state using a given rule
    with periodic boundary conditions.
    """
    n = len(state)
    # The rule is an 8-bit binary number.
    # It maps neighborhoods '111', '110', ..., '000' to the corresponding bits.
    rule_bin = format(rule_num, '08b')
    
    # We can create a mapping for clarity
    rule_map = {
        '111': rule_bin[0], '110': rule_bin[1], '101': rule_bin[2],
        '100': rule_bin[3], '011': rule_bin[4], '010': rule_bin[5],
        '001': rule_bin[6], '000': rule_bin[7]
    }
    
    next_state = ""
    for i in range(n):
        # Get neighbors with periodic (wrapping) boundary conditions
        left_neighbor = state[(i - 1 + n) % n]
        center_cell = state[i]
        right_neighbor = state[(i + 1) % n]
        
        neighborhood = left_neighbor + center_cell + right_neighbor
        
        # Append the new cell state based on the rule
        next_state += rule_map[neighborhood]
        
    return next_state

# Run the solver
solve_cellular_automaton()