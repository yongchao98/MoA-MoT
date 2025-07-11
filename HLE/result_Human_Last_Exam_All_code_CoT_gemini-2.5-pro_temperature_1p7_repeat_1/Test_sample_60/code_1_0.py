from collections import defaultdict

def parse_rules(rule_string):
    """
    Parses a compact Turing machine rule string into a nested dictionary.
    The format is {'State': {read_symbol: (next_state, write_symbol, move_direction)}}.
    """
    rules = {}
    states = ['A', 'B', 'C', 'D', 'E']
    for i in range(10):  # 10 rules total for 5 states and 2 symbols
        # Each rule is a 3-character chunk (e.g., 'B1L')
        rule_chunk = rule_string[i*3 : i*3 + 3]
        
        # Determine the current state and symbol this rule applies to
        state_name = states[i // 2]
        symbol_read = i % 2

        if state_name not in rules:
            rules[state_name] = {}

        # Extract the transition details
        next_state = rule_chunk[0]
        write_symbol = int(rule_chunk[1])
        move_direction = rule_chunk[2]

        rules[state_name][symbol_read] = (next_state, write_symbol, move_direction)
    return rules

def simulate_turing_machine(rule_string, machine_name, max_steps=500):
    """
    Simulates a Turing machine on an empty tape.
    Returns the number of steps until it halts or an indicator for non-halting.
    """
    rules = parse_rules(rule_string)
    tape = defaultdict(int)  # The tape is infinite and initialized with 0s
    state = 'A'              # Initial state
    head_pos = 0             # Initial head position
    steps = 0

    while state != 'H' and steps < max_steps:
        # Read the symbol on the tape at the current head position
        symbol_read = tape[head_pos]

        # Get the transition rule for the current state and symbol
        if state not in rules or symbol_read not in rules[state]:
            print(f"Machine {machine_name}: Undefined transition for state '{state}' reading '{symbol_read}'. Assuming non-halting.")
            return -1

        rule = rules[state][symbol_read]
        next_state, write_symbol, move = rule

        # Apply the rule: write to tape, change state, move head
        tape[head_pos] = write_symbol
        state = next_state
        if move == 'L':
            head_pos -= 1
        elif move == 'R':
            head_pos += 1
        
        steps += 1
    
    if steps >= max_steps:
        print(f"Machine {machine_name} did not halt within {max_steps} steps.")
        return -2  # Indicates a timeout
    
    print(f"Machine {machine_name} halted after {steps} steps.")
    return steps

def main():
    """
    Main function to run the simulations and find the winner.
    """
    # Define the rule strings for the three machines
    machines = [
        ("1", "B1LB0RC1RE0LA1LD0RC0RA1RC1LH1L"),
        ("2", "B1LC1RC1LB1LD1LE0RA1RD1RH1LA0R"),
        ("3", "B1LA1LC1RB1RA1LD1RA1LE1RH1LC0R")
    ]

    results = []
    print("--- Simulating Turing Machines ---")
    for name, rule_str in machines:
        steps = simulate_turing_machine(rule_str, name)
        # Store results only for machines that halt
        if steps > 0:
            results.append((name, steps))

    print("\n--- Final Analysis ---")
    if not results:
        print("None of the Turing Machines halted within the step limit.")
        return

    # Find the machine with the maximum number of steps among the ones that halted
    winner_name, max_steps = max(results, key=lambda item: item[1])

    print(f"Machine {winner_name} is the Turing Machine that halts after the most number of steps.")
    print(f"The number of steps is: {max_steps}")
    
    # Per instructions, output the final numerical result in the specified format
    print(f"\nFinal Result: Machine {winner_name} = {max_steps} steps")


if __name__ == "__main__":
    main()
<<<47>>>