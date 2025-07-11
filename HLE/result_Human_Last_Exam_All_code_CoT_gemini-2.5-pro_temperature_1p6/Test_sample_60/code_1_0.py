import collections

def simulate_turing_machine(rule_string):
    """
    Parses a rule string and simulates the Turing Machine.

    The rule string format is a space-separated list of transitions.
    For states A, B, C, D, E, the transitions for reading 0 and 1
    are listed in order.
    Example Rule: B1L -> Next State: B, Write: 1, Move: L
    """
    # 1. Parse the rules
    rules = {}
    parts = rule_string.split()
    states = ['A', 'B', 'C', 'D', 'E']
    part_index = 0
    for state in states:
        rules[state] = {}
        # Rule for symbol 0
        rule_part_0 = parts[part_index]
        rules[state][0] = (rule_part_0[0], int(rule_part_0[1]), rule_part_0[2])
        part_index += 1
        # Rule for symbol 1
        rule_part_1 = parts[part_index]
        rules[state][1] = (rule_part_1[0], int(rule_part_1[1]), rule_part_1[2])
        part_index += 1

    # 2. Initialize the simulation
    tape = collections.defaultdict(int)
    head = 0
    state = 'A'
    steps = 0
    max_steps = 100000  # Safety break for non-halting machines

    # 3. Run the simulation loop
    while state != 'H' and steps < max_steps:
        current_symbol = tape[head]
        
        # This machine definition does not have rules for state 'H'
        # so this check is sufficient.
        next_state, write_symbol, move = rules[state][current_symbol]

        tape[head] = write_symbol
        
        if move == 'L':
            head -= 1
        elif move == 'R':
            head += 1
        
        state = next_state
        steps += 1
        
    return steps if steps < max_steps else -1 # Return -1 if it didn't halt

def solve():
    """
    Solves the problem by simulating three Turing Machines and finding
    the one that halts after the most steps.
    """
    machines = [
        "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",
        "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",
        "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"
    ]
    
    results = []
    print("Simulating Turing Machines...")
    for i, m_str in enumerate(machines):
        num_steps = simulate_turing_machine(m_str)
        results.append(num_steps)
        print(f"Machine {i + 1} halts after {num_steps} steps.")

    max_steps = -1
    winning_machine_index = -1
    for i, steps in enumerate(results):
        if steps > max_steps:
            max_steps = steps
            winning_machine_index = i

    winning_machine_number = winning_machine_index + 1
    
    print("\n--- Final Answer ---")
    print(f"The Turing Machine that halts after the most number of steps is Machine {winning_machine_number}.")
    print(f"The number of steps is: {max_steps}")
    
    # Required final output format
    global answer_content
    answer_content = winning_machine_number


# Run the solver
solve()
# This will be captured by the platform. Do not copy it.
print(f'<<<{answer_content}>>>')