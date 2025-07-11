def parse_rules(rule_string):
    """
    Parses a compact rule string into a dictionary representing the TM's transition function.
    The format is assumed to be for states A, B, C, D, E and symbols 0, 1 in order.
    A rule "SWS'" means: Next State S, Write Symbol W, Move direction S'.
    """
    parts = rule_string.split()
    rules = {}
    states = "ABCDE"
    # The rule string contains 10 parts, for 5 states and 2 symbols (0, 1) each.
    # The order is (A,0), (A,1), (B,0), (B,1), ...
    for i, state_char in enumerate(states):
        rules[state_char] = {}
        # Rule for symbol 0
        rule_part_0 = parts[i * 2]
        rules[state_char][0] = {
            'next_state': rule_part_0[0],
            'write_symbol': int(rule_part_0[1]),
            'move': rule_part_0[2]
        }
        # Rule for symbol 1
        rule_part_1 = parts[i * 2 + 1]
        rules[state_char][1] = {
            'next_state': rule_part_1[0],
            'write_symbol': int(rule_part_1[1]),
            'move': rule_part_1[2]
        }
    return rules

def simulate_turing_machine(rules, max_steps=500000):
    """
    Simulates a Turing Machine given its rules, starting from an empty tape.
    Returns the number of steps to halt, or -1 if it doesn't halt within max_steps.
    """
    tape = {}  # A dictionary represents the sparse, infinite tape
    head_position = 0
    current_state = 'A'  # The start state is 'A'
    step_count = 0

    while step_count < max_steps:
        if current_state == 'H':
            return step_count  # Halt state reached

        # Read the symbol on the tape at the current head position. Default to 0 for blank cells.
        current_symbol = tape.get(head_position, 0)

        # Find the transition rule for the current state and symbol
        try:
            transition = rules[current_state][current_symbol]
        except KeyError:
            # This would happen if an unexpected state is reached (e.g., trying to execute from 'H')
            return -1

        # Apply the rule: write to tape, move head, and change state
        tape[head_position] = transition['write_symbol']
        
        if transition['move'] == 'R':
            head_position += 1
        elif transition['move'] == 'L':
            head_position -= 1
        
        current_state = transition['next_state']
        step_count += 1
    
    return -1 # Indicates non-halting within the step limit

# --- Main Execution ---

# The definitions for the three Turing Machines provided
tm_definitions = [
    "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",  # Machine 1
    "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",  # Machine 2
    "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"   # Machine 3
]

results = []
for i, definition in enumerate(tm_definitions):
    # Parse the rules for the current machine
    rules = parse_rules(definition)
    # Simulate the machine and store its step count
    steps = simulate_turing_machine(rules)
    results.append({'machine_number': i + 1, 'steps': steps})

# Filter out machines that did not halt
halting_machines = [res for res in results if res['steps'] >= 0]

if not halting_machines:
    print("None of the Turing Machines halted within the step limit.")
else:
    # Find the machine that took the most steps among those that halted
    winner = max(halting_machines, key=lambda item: item['steps'])
    
    # Print the final result including the machine number and step count
    print(f"Machine {winner['machine_number']} halts after the most number of steps, which is {winner['steps']}.")
