from collections import defaultdict

def simulate_turing_machine(name, rule_string):
    """
    Parses a rule string and simulates the corresponding Turing Machine.

    Args:
        name (str): The identifier for the machine (e.g., "Machine 1").
        rule_string (str): The string containing the machine's rules.

    Returns:
        int: The number of steps the machine takes to halt.
    """
    # 1. Parse the rule string into a usable dictionary structure
    rules = {}
    parts = rule_string.split()
    states = ['A', 'B', 'C', 'D', 'E']  # 5 states before Halt
    part_index = 0
    for state in states:
        rules[state] = {}
        # Rule for reading '0'
        rule_for_0 = parts[part_index]
        rules[state][0] = (rule_for_0[0], int(rule_for_0[1]), rule_for_0[2])
        part_index += 1
        # Rule for reading '1'
        rule_for_1 = parts[part_index]
        rules[state][1] = (rule_for_1[0], int(rule_for_1[1]), rule_for_1[2])
        part_index += 1

    # 2. Initialize the machine
    tape = defaultdict(int)  # Tape is initially all zeros
    head_position = 0
    current_state = 'A'
    steps = 0
    # Add a safety limit to prevent running forever on a non-halting machine
    max_steps = 1000000 

    # 3. Run the simulation loop
    while current_state != 'H' and steps < max_steps:
        steps += 1
        value_read = tape[head_position]
        
        # Get the rule for the current state and tape value
        next_state, value_to_write, move_direction = rules[current_state][value_read]
        
        # Update tape, state, and head position
        tape[head_position] = value_to_write
        current_state = next_state
        
        if move_direction == 'R':
            head_position += 1
        elif move_direction == 'L':
            head_position -= 1

    # 4. Return the total steps
    return steps

def main():
    """
    Main function to run simulations and find the winner.
    """
    machines = {
        1: "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",
        2: "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",
        3: "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"
    }

    results = {}
    print("Simulating Turing Machines...")
    print("-" * 30)

    for i, rules in machines.items():
        num_steps = simulate_turing_machine(f"Machine {i}", rules)
        results[i] = num_steps
        # Per the instructions, output the step count for each machine
        print(f"Machine {i} halts after {num_steps} steps.")

    print("-" * 30)
    
    # Find the machine with the most steps
    if not results:
        print("No results to analyze.")
        return

    winner_machine_num = max(results, key=results.get)
    max_steps = results[winner_machine_num]

    print(f"The Turing Machine that halts after the most number of steps is Machine {winner_machine_num}.")
    print(f"The number of steps is {max_steps}.")

if __name__ == "__main__":
    main()
<<<2>>>