import sys

def solve_turing_machine_race():
    """
    Parses and simulates three Turing machines to find the one that halts
    after the most steps, starting from an empty tape.
    """

    # The rules for the three Turing Machines provided by the user.
    machine_rules_str = [
        "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",  # Machine 1
        "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",  # Machine 2
        "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"   # Machine 3
    ]

    def parse_rules(rule_string):
        """Parses a flat string of rules into a nested dictionary."""
        parts = rule_string.split()
        states = "ABCDE"
        rules = {}
        rule_index = 0
        for state in states:
            rules[state] = {}
            for symbol in [0, 1]:
                rule_spec = parts[rule_index]
                new_state = rule_spec[0]
                new_symbol = int(rule_spec[1])
                # Convert 'L' to -1 and 'R' to 1 for head movement
                direction = -1 if rule_spec[2] == 'L' else 1
                rules[state][symbol] = (new_state, new_symbol, direction)
                rule_index += 1
        return rules

    def simulate_tm(rules, max_steps=500000):
        """Simulates a Turing Machine and returns the number of steps to halt."""
        tape = {}
        head_pos = 0
        state = 'A'
        steps = 0

        while state != 'H':
            if steps >= max_steps:
                return -1  # Indicates non-halting or exceeding step limit

            current_symbol = tape.get(head_pos, 0)
            
            try:
                new_state, new_symbol, direction = rules[state][current_symbol]
            except KeyError:
                # This handles cases where a state is referenced that isn't defined.
                # In this problem, we assume only 'H' is a valid halting state not in the ruleset.
                print(f"Error: Machine halted due to undefined rule for State {state}, Symbol {current_symbol}", file=sys.stderr)
                return steps 

            tape[head_pos] = new_symbol
            state = new_state
            head_pos += direction
            steps += 1
            
        return steps

    # --- Main execution logic ---
    max_steps = -1
    winner_index = -1
    
    print("--- Turing Machine Simulation ---")
    for i, rule_string in enumerate(machine_rules_str):
        parsed_rules = parse_rules(rule_string)
        steps = simulate_tm(parsed_rules)
        machine_num = i + 1
        
        if steps == -1:
            print(f"Machine {machine_num}: Did not halt within the step limit.")
        else:
            print(f"Machine {machine_num}: Halted in {steps} steps.")
            if steps > max_steps:
                max_steps = steps
                winner_index = i
    
    print("\n--- Result ---")
    if winner_index == -1:
        print("No machine that was tested halted within the step limit.")
    else:
        winner_name = winner_index + 1
        print(f"Machine {winner_name} halts after the most number of steps.")
        print(f"The number of steps is: {max_steps}")

# Run the simulation
solve_turing_machine_race()