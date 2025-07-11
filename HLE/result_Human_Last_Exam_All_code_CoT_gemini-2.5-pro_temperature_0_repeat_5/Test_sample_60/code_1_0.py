import collections

def solve_turing_machines():
    """
    Parses, simulates, and compares three Turing Machines to find the one
    that halts after the most steps on an empty tape.
    """

    tm_definitions = [
        "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",  # Machine 1
        "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",  # Machine 2
        "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"   # Machine 3
    ]

    def parse_rules(definition_string):
        """Parses a rule string into a dictionary."""
        rules = {}
        states = ['A', 'B', 'C', 'D', 'E']
        parts = definition_string.split()
        i = 0
        for state in states:
            for read_val in [0, 1]:
                rule_part = parts[i]
                new_state = rule_part[0]
                write_val = int(rule_part[1])
                move = rule_part[2]
                rules[(state, read_val)] = (write_val, move, new_state)
                i += 1
        return rules

    def simulate(rules):
        """Simulates a Turing Machine and returns the number of steps to halt."""
        tape = collections.defaultdict(int)
        position = 0
        state = 'A'
        steps = 0
        # Set a high limit to catch potential infinite loops
        max_steps_limit = 100000

        while state != 'H':
            if steps >= max_steps_limit:
                return f"Did not halt within {max_steps_limit} steps"

            read_val = tape[position]
            
            if (state, read_val) not in rules:
                # Implicit halt if no rule is defined, though the problem implies 'H' is the only halt state.
                break

            write_val, move, new_state = rules[(state, read_val)]
            
            tape[position] = write_val
            
            if move == 'R':
                position += 1
            elif move == 'L':
                position -= 1
            
            state = new_state
            steps += 1
            
        return steps

    results = []
    for i, definition in enumerate(tm_definitions):
        rules = parse_rules(definition)
        steps = simulate(rules)
        results.append(steps)
        print(f"Machine {i + 1} halts after {steps} steps.")

    max_steps = -1
    winner_index = -1
    for i, steps in enumerate(results):
        if isinstance(steps, int) and steps > max_steps:
            max_steps = steps
            winner_index = i

    if winner_index != -1:
        print(f"\nMachine {winner_index + 1} halts after the most number of steps, which is {max_steps}.")
    else:
        print("\nCould not determine a winner.")

solve_turing_machines()