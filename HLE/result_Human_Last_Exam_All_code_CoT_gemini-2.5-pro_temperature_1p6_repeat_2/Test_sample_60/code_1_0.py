import collections

def solve_turing_machine_race():
    """
    Simulates three Turing Machines to find which one halts after the most steps.
    """

    machines = [
        "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",
        "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",
        "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"
    ]

    def simulate(rule_string):
        """
        Simulates a single Turing Machine based on its rule string.
        """
        # --- 1. Parse the rules ---
        rules = {}
        parts = rule_string.split()
        states_in_order = "ABCDE"
        for i, state_char in enumerate(states_in_order):
            # Rule for reading '0'
            rule_for_0 = parts[i * 2]
            rules[(state_char, 0)] = (rule_for_0[0], int(rule_for_0[1]), rule_for_0[2])
            # Rule for reading '1'
            rule_for_1 = parts[i * 2 + 1]
            rules[(state_char, 1)] = (rule_for_1[0], int(rule_for_1[1]), rule_for_1[2])

        # --- 2. Initialize the machine ---
        tape = collections.defaultdict(int)
        head_position = 0
        current_state = 'A'
        steps = 0
        
        # --- 3. Run the simulation ---
        # Added a safety break for extremely long or non-halting machines
        max_steps = 100000 
        while current_state != 'H' and steps < max_steps:
            steps += 1
            symbol_read = tape[head_position]
            
            new_state, symbol_to_write, move_direction = rules[(current_state, symbol_read)]
            
            tape[head_position] = symbol_to_write
            current_state = new_state
            
            if move_direction == 'R':
                head_position += 1
            else:  # 'L'
                head_position -= 1

        return steps

    # --- 4. Find the machine with the most steps ---
    max_steps = -1
    winner_index = -1
    
    for i, machine_rules in enumerate(machines):
        num_steps = simulate(machine_rules)
        if num_steps > max_steps:
            max_steps = num_steps
            winner_index = i + 1

    # --- 5. Print the result ---
    print(f"Machine number {winner_index} halts after the most steps.")
    print(f"The number of steps is {max_steps}.")

solve_turing_machine_race()
<<<Machine number 3 halts after the most steps.
The number of steps is 96.>>>