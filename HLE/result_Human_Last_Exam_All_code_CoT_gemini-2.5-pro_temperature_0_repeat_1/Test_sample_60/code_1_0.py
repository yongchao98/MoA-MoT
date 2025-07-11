import collections

def simulate_turing_machine(rule_string: str, max_steps: int = 100000):
    """
    Simulates a Turing machine based on a rule string.

    The rule string format is a concatenation of rules for states A, B, C, ...
    Each state has two rules (for reading 0 and 1), each 3 characters long.
    Example rule: 'B1L' means "move to state B, write 1, move Left".

    Args:
        rule_string: The string defining the machine's rules.
        max_steps: A safeguard limit to prevent infinite loops.

    Returns:
        The number of steps until the machine halts, or -1 if it exceeds max_steps.
    """
    # Parse the rule string into a more usable dictionary format
    rules = {}
    state_chars = "ABCDEFGH"
    # Remove spaces for easier parsing
    cleaned_rule_string = rule_string.replace(" ", "")
    
    # Each state (A, B, C...) has two 3-character rules, making a 6-character block
    rule_blocks = [cleaned_rule_string[i:i+6] for i in range(0, len(cleaned_rule_string), 6)]

    for i, block in enumerate(rule_blocks):
        current_state = state_chars[i]
        rules[current_state] = {}
        
        # Rule for reading a 0
        rule_for_0 = block[0:3]
        next_state_0, write_val_0, move_0 = rule_for_0[0], int(rule_for_0[1]), rule_for_0[2]
        rules[current_state][0] = (next_state_0, write_val_0, 1 if move_0 == 'R' else -1)
        
        # Rule for reading a 1
        rule_for_1 = block[3:6]
        next_state_1, write_val_1, move_1 = rule_for_1[0], int(rule_for_1[1]), rule_for_1[2]
        rules[current_state][1] = (next_state_1, write_val_1, 1 if move_1 == 'R' else -1)

    # Initialize the simulation
    tape = collections.defaultdict(int)  # Tape is infinite and defaults to 0
    state = 'A'
    head_position = 0
    steps = 0

    # Run the simulation loop
    while state != 'H':
        if steps >= max_steps:
            return -1  # Machine did not halt within the step limit

        current_tape_value = tape[head_position]
        
        # Get the rule for the current state and tape value
        next_state, value_to_write, move_direction = rules[state][current_tape_value]
        
        # Apply the rule
        tape[head_position] = value_to_write
        head_position += move_direction
        state = next_state
        steps += 1
        
    return steps

def main():
    """
    Main function to define, run, and compare the Turing machines.
    """
    machines = [
        "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",
        "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",
        "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"
    ]

    results = []
    for i, machine_code in enumerate(machines):
        num_steps = simulate_turing_machine(machine_code)
        results.append((i + 1, num_steps))
        if num_steps != -1:
            print(f"Machine {i + 1} halts after {num_steps} steps.")
        else:
            print(f"Machine {i + 1} did not halt within the step limit.")

    # Find the machine with the most steps
    if not results:
        print("\nNo simulations were run.")
        return

    # Filter out non-halting machines for comparison
    halting_machines = [res for res in results if res[1] != -1]

    if not halting_machines:
        print("\nNone of the machines halted within the step limit.")
        return

    winner = max(halting_machines, key=lambda item: item[1])
    winner_machine_number = winner[0]
    winner_steps = winner[1]

    print(f"\nMachine {winner_machine_number} halts after the most number of steps.")
    print(f"The number of steps is: {winner_steps}")
    
    # The final answer format as requested
    print(f"\n<<<Machine {winner_machine_number}>>>")


if __name__ == "__main__":
    main()