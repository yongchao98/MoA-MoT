import collections

def simulate_turing_machine(rules_string):
    """
    Simulates a Turing Machine based on a given rule string.
    
    The machine starts in state 'A' on an empty tape (all zeros) at position 0.
    The simulation runs until it reaches the halt state 'H'.
    
    Args:
        rules_string: A space-separated string of rules.
    
    Returns:
        The total number of steps taken to halt. Returns -1 if it exceeds a max step limit.
    """
    # 1. Parse the rules from the input string
    rules = {}
    states = ['A', 'B', 'C', 'D', 'E']
    rule_parts = rules_string.split()
    
    rule_index = 0
    for state in states:
        # Rule for reading symbol 0
        rule0_str = rule_parts[rule_index]
        rule0 = (rule0_str[0], int(rule0_str[1]), rule0_str[2]) # (next_state, write_symbol, move_direction)
        
        # Rule for reading symbol 1
        rule1_str = rule_parts[rule_index + 1]
        rule1 = (rule1_str[0], int(rule1_str[1]), rule1_str[2])
        
        rules[state] = [rule0, rule1]
        rule_index += 2

    # 2. Initialize the machine's state
    tape = collections.defaultdict(int)  # The tape is infinite and defaults to 0
    head_position = 0
    current_state = 'A'
    steps = 0
    
    # Safety limit to prevent hanging on a non-halting machine
    max_steps = 200000 

    # 3. Run the simulation loop
    while current_state != 'H':
        if steps >= max_steps:
             print(f"Warning: Machine exceeded {max_steps} steps and is assumed to be non-halting.")
             return -1

        current_symbol = tape[head_position]
        
        # Get the rule for the current state and symbol
        next_state, write_symbol, move_direction = rules[current_state][current_symbol]
        
        # Update the tape
        tape[head_position] = write_symbol
        
        # Update the state
        current_state = next_state
        
        # Move the head
        if move_direction == 'R':
            head_position += 1
        elif move_direction == 'L':
            head_position -= 1
            
        steps += 1
        
    return steps

def main():
    """
    Main function to run simulations for all machines and find the one with the most steps.
    """
    turing_machines = [
        "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",  # Machine 1
        "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",  # Machine 2
        "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"   # Machine 3
    ]

    results = []
    print("--- Simulating Turing Machines ---")
    for i, definition in enumerate(turing_machines):
        machine_id = i + 1
        num_steps = simulate_turing_machine(definition)
        results.append(num_steps)
        print(f"Machine {machine_id} halted after {num_steps} steps.")

    max_steps = -1
    winner_id = -1
    for i, steps in enumerate(results):
        if steps > max_steps:
            max_steps = steps
            winner_id = i + 1
            
    print("\n--- Final Result ---")
    if winner_id != -1:
        print(f"The Turing Machine that halts after the most number of steps is Machine {winner_id}.")
        print(f"It halts after {max_steps} steps.")
    else:
        print("Could not determine a winner.")

if __name__ == "__main__":
    main()