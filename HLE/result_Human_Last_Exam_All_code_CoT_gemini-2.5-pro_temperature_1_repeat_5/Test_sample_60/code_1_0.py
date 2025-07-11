import collections

def simulate_turing_machine(rule_string):
    """
    Simulates a 5-state, 2-symbol Turing machine based on a rule string.

    The machine starts in state 'A' with the head at position 0 on an
    infinite tape of 0s. The simulation halts when state 'H' is reached.
    """
    # Parse the rule string into a more usable format (a dictionary)
    parts = rule_string.split()
    states_order = ['A', 'B', 'C', 'D', 'E']
    rules = collections.defaultdict(dict)
    for i, state_char in enumerate(states_order):
        # Rule for reading 0
        rule_0_str = parts[i * 2]
        rules[state_char][0] = (rule_0_str[0], int(rule_0_str[1]), rule_0_str[2])
        # Rule for reading 1
        rule_1_str = parts[i * 2 + 1]
        rules[state_char][1] = (rule_1_str[0], int(rule_1_str[1]), rule_1_str[2])

    # --- Simulation Setup ---
    # Use a defaultdict for a tape that is infinite in both directions
    # and defaults to 0 for any unvisited cell.
    tape = collections.defaultdict(int)
    head_position = 0
    current_state = 'A'  # Start state
    steps = 0
    
    # A safety limit to prevent running forever on non-halting machines.
    max_steps = 50000 

    while current_state != 'H' and steps < max_steps:
        # Read the symbol at the current head position
        symbol_read = tape[head_position]

        # Find the appropriate rule for the current state and symbol
        rule = rules[current_state][symbol_read]
        new_state, symbol_to_write, direction = rule

        # Apply the rule:
        # 1. Write the new symbol to the tape
        tape[head_position] = symbol_to_write

        # 2. Move the head
        if direction == 'L':
            head_position -= 1
        elif direction == 'R':
            head_position += 1

        # 3. Update the state
        current_state = new_state
        
        # 4. Count the step
        steps += 1
        
    if current_state == 'H':
        return steps
    else:
        # Return -1 if the machine did not halt within the step limit
        return -1

# --- Main Execution ---

# Define the three Turing Machines from the problem description
machines = {
    "Machine 1": "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",
    "Machine 2": "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",
    "Machine 3": "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"
}

# Run the simulation for each machine and store the results
results = {}
for name, rule_string in machines.items():
    steps = simulate_turing_machine(rule_string)
    results[name] = steps

# Find the machine that ran for the most steps
winner_name = max(results, key=results.get)
max_steps = results[winner_name]

# Print the results for all machines and announce the winner
print("--- Turing Machine Simulation Results ---")
for name, step_count in results.items():
    if step_count != -1:
        print(f"{name} halted after {step_count} steps.")
    else:
        print(f"{name} did not halt within the step limit.")

print("\n--- Conclusion ---")
print(f"The Turing Machine that halts after the most number of steps is {winner_name}.")
print(f"The number of steps is: {max_steps}")