import collections

def simulate_turing_machine(rules_string):
    """
    Simulates a Turing Machine based on a given rule string.

    Args:
        rules_string: A string containing the machine's transition rules.

    Returns:
        The number of steps the machine took to halt, or -1 if it doesn't halt
        within the maximum step limit.
    """
    states = ['A', 'B', 'C', 'D', 'E']
    rules_map = {}
    rule_parts = rules_string.split()
    
    # Parse the compact rule string into a nested dictionary for easy lookup
    part_idx = 0
    for state_char in states:
        rules_map[state_char] = {}
        for read_val in [0, 1]:
            rule_str = rule_parts[part_idx]
            new_state = rule_str[0]
            write_val = int(rule_str[1])
            move = rule_str[2]
            rules_map[state_char][read_val] = (new_state, write_val, move)
            part_idx += 1

    # Initialize the machine's components
    # A defaultdict is used for the tape, so unwritten cells default to 0.
    tape = collections.defaultdict(int)
    head = 0
    state = 'A'  # Start state is always 'A'
    steps = 0
    # A safety limit to prevent true infinite loops from crashing the script
    max_steps = 10000

    # Main simulation loop
    while state != 'H':
        if steps >= max_steps:
            return -1  # Assume it does not halt in a reasonable time

        # Read the value from the tape at the current head position
        current_val_on_tape = tape[head]

        # Get the rule for the current state and tape value
        new_state, write_val, move_dir = rules_map[state][current_val_on_tape]

        # Apply the rule
        tape[head] = write_val
        if move_dir == 'R':
            head += 1
        else:  # 'L'
            head -= 1
        
        state = new_state
        steps += 1

    return steps

def main():
    """
    Main function to run the simulations and find the winner.
    """
    machine_rules = [
        "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",
        "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",
        "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R",
    ]

    results = []
    for rules in machine_rules:
        num_steps = simulate_turing_machine(rules)
        results.append(num_steps)

    max_steps = -1
    winner_index = -1
    for i, steps in enumerate(results):
        if steps > max_steps:
            max_steps = steps
            winner_index = i

    print(f"Machine 1 halted after {results[0]} steps.")
    print(f"Machine 2 halted after {results[1]} steps.")
    print(f"Machine 3 halted after {results[2]} steps.")
    
    if winner_index != -1:
        print(f"\nThe Turing Machine that halts after the most number of steps is machine {winner_index + 1}.")
        print(f"It halts after {max_steps} steps.")
    else:
        print("\nCould not determine a winner.")

if __name__ == "__main__":
    main()