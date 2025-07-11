import collections

def solve_turing_machine_problem():
    """
    Parses, simulates, and finds the Turing Machine that halts after the most steps.
    """

    tm_definitions = [
        "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",  # Machine 1
        "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",  # Machine 2
        "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R",  # Machine 3
    ]

    def parse_rules(rule_string):
        """Converts a compact rule string into a dictionary."""
        parts = rule_string.split()
        states_order = ['A', 'B', 'C', 'D', 'E']
        rules = {}
        part_index = 0
        for state in states_order:
            rules[state] = {}
            # Rule for reading 0
            rule_part_0 = parts[part_index]
            rules[state][0] = (rule_part_0[0], int(rule_part_0[1]), rule_part_0[2])
            part_index += 1
            # Rule for reading 1
            rule_part_1 = parts[part_index]
            rules[state][1] = (rule_part_1[0], int(rule_part_1[1]), rule_part_1[2])
            part_index += 1
        return rules

    def simulate_tm(rules):
        """Simulates a Turing Machine and returns the number of steps to halt."""
        # A defaultdict simplifies tape management, assuming 0 for unwritten cells.
        tape = collections.defaultdict(int)
        head_position = 0
        current_state = 'A'
        steps = 0
        
        # Safety break to prevent running forever in case of a non-halting machine.
        max_steps = 1000000 

        while current_state != 'H':
            if steps >= max_steps:
                print("Warning: Exceeded maximum steps, assuming non-halting.")
                return -1

            read_value = tape[head_position]
            
            # Get the transition rule for the current state and tape value.
            next_state, write_value, move_direction = rules[current_state][read_value]
            
            # Update the tape.
            tape[head_position] = write_value
            
            # Move the head.
            if move_direction == 'R':
                head_position += 1
            elif move_direction == 'L':
                head_position -= 1
            
            # Transition to the next state.
            current_state = next_state
            
            # Increment the step counter.
            steps += 1
            
        return steps

    step_counts = []
    for definition in tm_definitions:
        parsed_rules = parse_rules(definition)
        steps = simulate_tm(parsed_rules)
        step_counts.append(steps)
    
    # Find the machine with the maximum number of steps.
    max_steps = -1
    winner_index = -1
    if step_counts:
        max_steps = max(step_counts)
        winner_index = step_counts.index(max_steps)

    # Print the final result including the winning machine number and its step count.
    winning_machine_number = winner_index + 1
    print(f"Machine {winning_machine_number} halts after the most number of steps.")
    print(f"Number of steps: {max_steps}")

solve_turing_machine_problem()