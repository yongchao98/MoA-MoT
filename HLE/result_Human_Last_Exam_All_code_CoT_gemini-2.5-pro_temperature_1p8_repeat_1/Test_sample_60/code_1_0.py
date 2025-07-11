import collections

def solve():
    """
    Solves the Turing Machine halting problem for the given machines by simulating them.
    """

    def parse_rules(rule_string):
        """Parses the compact rule string into a dictionary for the simulator."""
        rules = {}
        parts = rule_string.split()
        states = "ABCDE"
        for i, state in enumerate(states):
            # Rule for reading 0
            rule0_str = parts[i * 2]
            new_state0 = rule0_str[0]
            write_val0 = int(rule0_str[1])
            move_dir0 = rule0_str[2]
            rules[(state, 0)] = (new_state0, write_val0, move_dir0)

            # Rule for reading 1
            rule1_str = parts[i * 2 + 1]
            new_state1 = rule1_str[0]
            write_val1 = int(rule1_str[1])
            move_dir1 = rule1_str[2]
            rules[(state, 1)] = (new_state1, write_val1, move_dir1)
        return rules

    def simulate_turing_machine(rule_string, max_steps=1000000):
        """
        Simulates a Turing Machine based on a rule string.
        Returns the number of steps until it halts.
        """
        rules = parse_rules(rule_string)
        # defaultdict(int) conveniently models an infinite tape of 0s
        tape = collections.defaultdict(int)
        head_pos = 0
        state = 'A'  # All machines start at state A
        steps = 0

        while state != 'H':
            if steps >= max_steps:
                return f"Did not halt within {max_steps} steps"

            current_val = tape[head_pos]
            
            # Look up rule, unpack new state, write value, and move direction
            new_state, write_val, move_dir = rules[(state, current_val)]

            # Apply the rule
            tape[head_pos] = write_val

            if move_dir == 'L':
                head_pos -= 1
            elif move_dir == 'R':
                head_pos += 1

            state = new_state
            steps += 1
            
        return steps

    tm_definitions = {
        1: "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",
        2: "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",
        3: "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"
    }

    results = {}
    print("Simulating Turing Machines...")
    for machine_id, rule_str in tm_definitions.items():
        step_count = simulate_turing_machine(rule_str)
        results[machine_id] = step_count
        print(f"Machine {machine_id} halted after {step_count} steps.")

    # Find the machine with the most steps from the valid (integer) results
    valid_results = {k: v for k, v in results.items() if isinstance(v, int)}
    
    if not valid_results:
        print("No machine halted within the step limit.")
        return
    
    # Determine the winner
    winner_id = max(valid_results, key=valid_results.get)
    max_steps = valid_results[winner_id]
    
    print("\n--- Final Result ---")
    print(f"The Turing Machine that halts after the most steps is Machine {winner_id}.")
    print(f"It halts after {max_steps} steps.")

solve()
<<<3>>>