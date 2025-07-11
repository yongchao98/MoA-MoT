import sys
from collections import defaultdict

def parse_rules(rule_string):
    """
    Parses a compact Turing Machine rule string into a usable dictionary.
    Format assumes states 'A' through 'E' and rules for reading 0 and 1.
    e.g., "B1L B0R..." means State A, read 0 -> (B, 1, L); read 1 -> (B, 0, R)
    """
    states = ['A', 'B', 'C', 'D', 'E']
    rule_parts = rule_string.split()
    
    rules = {}
    if len(rule_parts) != len(states) * 2:
        raise ValueError("Invalid rule string format")

    for i, state_char in enumerate(states):
        # Rule for reading '0' is at index 2*i
        rule_for_0_str = rule_parts[i * 2]
        # Rule for reading '1' is at index 2*i + 1
        rule_for_1_str = rule_parts[i * 2 + 1]
        
        # Parse the rule for reading '0'
        next_state_0 = rule_for_0_str[0]
        write_symbol_0 = int(rule_for_0_str[1])
        move_dir_0 = rule_for_0_str[2]
        
        # Parse the rule for reading '1'
        next_state_1 = rule_for_1_str[0]
        write_symbol_1 = int(rule_for_1_str[1])
        move_dir_1 = rule_for_1_str[2]
        
        rules[state_char] = {
            0: (next_state_0, write_symbol_0, move_dir_0),
            1: (next_state_1, write_symbol_1, move_dir_1)
        }
    return rules

def simulate_turing_machine(rules, max_steps=200000):
    """
    Simulates a Turing Machine given a set of rules.
    """
    tape = defaultdict(int)  # The tape is infinite and defaults to 0
    head_position = 0
    current_state = 'A'
    steps = 0

    while current_state != 'H' and steps < max_steps:
        # Read the symbol from the tape
        read_symbol = tape[head_position]

        # Get the appropriate rule
        try:
            next_state, write_symbol, move_dir = rules[current_state][read_symbol]
        except KeyError:
            print(f"Error: No rule found for state '{current_state}' reading '{read_symbol}'. Halting.")
            return -1

        # Execute the rule
        tape[head_position] = write_symbol

        if move_dir == 'R':
            head_position += 1
        elif move_dir == 'L':
            head_position -= 1

        current_state = next_state
        steps += 1

    if current_state == 'H':
        return steps
    else:
        # Machine did not halt within the maximum allowed steps
        return -1

def main():
    """
    Main function to parse, simulate, and compare the Turing Machines.
    """
    tm_definitions = [
        "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",
        "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",
        "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"
    ]

    step_counts = []
    for i, tm_def in enumerate(tm_definitions):
        try:
            rules = parse_rules(tm_def)
            steps = simulate_turing_machine(rules)
            step_counts.append(steps)
            if steps == -1:
                print(f"Machine {i + 1} did not halt within the step limit.")
            else:
                print(f"Machine {i + 1} halts after {steps} steps.")
        except (ValueError, IndexError) as e:
            print(f"Error processing Machine {i + 1}: {e}")
            step_counts.append(-1)
            
    if not any(s > 0 for s in step_counts):
        print("None of the machines halted.")
        return

    max_steps = -1
    winning_machine_index = -1
    for i, steps in enumerate(step_counts):
        if steps > max_steps:
            max_steps = steps
            winning_machine_index = i

    print("\n--- Final Result ---")
    print(f"The Turing Machine that halts after the most number of steps is Machine {winning_machine_index + 1}.")
    print(f"The number of steps is: {max_steps}")

if __name__ == "__main__":
    main()