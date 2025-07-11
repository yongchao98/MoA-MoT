import sys

def apply_rule(state, rule_num):
    """
    Applies an elementary cellular automaton rule to a state, assuming
    periodic boundary conditions.
    """
    # An elementary rule is defined by its 8-bit binary representation.
    # The bits correspond to the output for neighborhoods 111, 110, ..., 000.
    try:
        rule_bin = format(rule_num, '08b')
    except (ValueError, TypeError):
        return None # Should not happen with integer rule_num

    rule_map = {
        '111': rule_bin[0], '110': rule_bin[1], '101': rule_bin[2],
        '100': rule_bin[3], '011': rule_bin[4], '010': rule_bin[5],
        '001': rule_bin[6], '000': rule_bin[7]
    }
    
    n = len(state)
    next_state = []
    
    # Iterate through each cell to calculate its next state.
    for i in range(n):
        # Periodic boundaries mean the grid wraps around.
        # The left neighbor of the first cell is the last cell.
        # The right neighbor of the last cell is the first cell.
        left_neighbor = state[i - 1]
        current_cell = state[i]
        right_neighbor = state[(i + 1) % n]
        
        neighborhood = left_neighbor + current_cell + right_neighbor
        
        next_state.append(rule_map[neighborhood])
        
    return "".join(next_state)

def find_intermediate_step():
    """
    Searches through all 256 elementary CA rules to find the unique
    intermediate step between the two given generations.
    """
    initial_state = "01101001"
    final_state = "10000111"
    
    valid_solutions = []

    # Iterate through all 256 possible elementary rules.
    for rule in range(256):
        # Step 1: Calculate the potential intermediate state.
        intermediate_state = apply_rule(initial_state, rule)
        
        # Step 2: Calculate the next state from the intermediate one.
        calculated_final_state = apply_rule(intermediate_state, rule)
        
        # Check if the result matches the given final state.
        if calculated_final_state == final_state:
            valid_solutions.append(intermediate_state)
            
    # Use a set to find the unique intermediate states found.
    unique_solutions = set(valid_solutions)
    
    if len(unique_solutions) == 1:
        # The problem statement implies there is only one valid solution.
        print(unique_solutions.pop())
    elif len(unique_solutions) == 0:
        print("No valid solution found under standard periodic boundary conditions.", file=sys.stderr)
    else:
        # This case is not expected based on the problem statement.
        print(f"Error: Found {len(unique_solutions)} unique solutions.", file=sys.stderr)

find_intermediate_step()