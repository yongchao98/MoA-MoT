import collections

def parse_rules(rule_string):
    """Parses the compact rule string into a usable dictionary."""
    parts = rule_string.split()
    rules = collections.defaultdict(dict)
    # States are A, B, C, D, E in alphabetical order.
    # H is the Halt state.
    states = ['A', 'B', 'C', 'D', 'E']
    rule_idx = 0
    for state in states:
        # Rule for reading symbol 0
        rule_part_0 = parts[rule_idx]
        rules[state][0] = (rule_part_0[0], int(rule_part_0[1]), rule_part_0[2])
        rule_idx += 1
        # Rule for reading symbol 1
        rule_part_1 = parts[rule_idx]
        rules[state][1] = (rule_part_1[0], int(rule_part_1[1]), rule_part_1[2])
        rule_idx += 1
    return rules

def simulate_tm(rules, max_steps=500000):
    """Simulates a Turing Machine based on the given rules."""
    tape = collections.defaultdict(int)
    head_pos = 0
    state = 'A' # Initial state
    steps = 0

    while state != 'H' and steps < max_steps:
        current_symbol = tape[head_pos]
        
        # Get the rule for the current state and symbol
        new_state, write_symbol, move_dir = rules[state][current_symbol]

        # Apply the rule
        tape[head_pos] = write_symbol

        if move_dir == 'L':
            head_pos -= 1
        elif move_dir == 'R':
            head_pos += 1

        state = new_state
        steps += 1
    
    if state == 'H':
        return "Halts", steps
    else:
        return f"Does not halt within {max_steps} steps", steps

def find_longest_halting_machine():
    """
    Parses, simulates, and finds the Turing Machine that halts after the most steps.
    """
    tm_definitions = [
        "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",
        "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",
        "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"
    ]

    results = []
    print("Simulating Turing Machines...")
    for i, tm_def in enumerate(tm_definitions, 1):
        rules = parse_rules(tm_def)
        status, steps = simulate_tm(rules)
        results.append({'machine': i, 'steps': steps, 'status': status})
        print(f"Machine {i}: {status} in {steps} steps.")

    max_steps = -1
    best_machine_info = None

    for res in results:
        # We only care about machines that actually halt
        if res['status'] == 'Halts' and res['steps'] > max_steps:
            max_steps = res['steps']
            best_machine_info = res

    print("\n--- Analysis Complete ---")
    if best_machine_info is None:
        print("None of the Turing Machines halted within the step limit.")
    else:
        best_machine_number = best_machine_info['machine']
        step_count = best_machine_info['steps']
        print(f"The machine that halts after the most steps is Machine {best_machine_number}.")
        # Final result output as requested
        print(f"Machine {best_machine_number} halts after {step_count} steps.")

if __name__ == '__main__':
    find_longest_halting_machine()