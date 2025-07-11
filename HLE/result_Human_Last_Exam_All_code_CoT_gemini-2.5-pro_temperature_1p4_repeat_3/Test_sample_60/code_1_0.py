def simulate_turing_machine(rules_str):
    """
    Simulates a Turing Machine based on a given rule string.

    Args:
        rules_str: A string containing the machine's transition rules.

    Returns:
        The total number of steps until the machine halts.
    """
    # 1. Parse the rules from the input string
    rules = {}
    states = "ABCDE"
    parts = rules_str.split()
    for i, state_char in enumerate(states):
        rules[state_char] = {}
        # Rule for reading '0'
        rule_for_0 = parts[i * 2]
        rules[state_char][0] = (rule_for_0[0], int(rule_for_0[1]), rule_for_0[2])
        # Rule for reading '1'
        rule_for_1 = parts[i * 2 + 1]
        rules[state_char][1] = (rule_for_1[0], int(rule_for_1[1]), rule_for_1[2])

    # 2. Initialize the machine
    tape = {}  # Using a dictionary to represent the sparse infinite tape
    head_pos = 0
    current_state = 'A'
    steps = 0
    # A safeguard against non-halting machines
    max_steps = 100000

    # 3. Run the simulation loop
    while current_state != 'H' and steps < max_steps:
        steps += 1
        # Read the value from the tape at the current head position.
        # .get(head_pos, 0) defaults to 0 if the cell hasn't been written to.
        current_val = tape.get(head_pos, 0)

        # Get the rule for the current state and tape value
        next_state, write_val, direction = rules[current_state][current_val]

        # Update tape, state, and head position based on the rule
        tape[head_pos] = write_val
        current_state = next_state
        if direction == 'R':
            head_pos += 1
        elif direction == 'L':
            head_pos -= 1

    return steps

if __name__ == '__main__':
    # Define the rule strings for the three machines
    machines = {
        1: "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",
        2: "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",
        3: "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"
    }

    results = {}
    print("Simulating Turing Machines...")
    for i, tm_str in machines.items():
        num_steps = simulate_turing_machine(tm_str)
        results[i] = num_steps
        # Outputting the number of steps for each machine
        print(f"Machine {i} halts after {num_steps} steps.")

    # Find the machine that ran for the most steps
    winner_machine = max(results, key=results.get)
    max_steps = results[winner_machine]

    print("\n--- Result ---")
    print(f"Machine {winner_machine} halts after the most number of steps.")
    print(f"The number of steps is: {max_steps}")