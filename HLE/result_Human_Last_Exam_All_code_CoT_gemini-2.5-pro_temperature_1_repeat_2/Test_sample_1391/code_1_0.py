import sys

def simulate_turing_machine():
    """
    Simulates a Turing machine based on a predefined program and initial tape state.
    """
    # 1. Initialization
    # Define the Turing machine's program
    instructions_str = """
0 1 0 r 0
0 0 0 r 0
0 _ _ l 1
1 0 0 l 1
1 1 1 l 1
1 _ _ r halt
""".strip().split('\n')

    # Define the starting conditions
    start_tape_str = "1H10"
    start_state = '0'

    # Parse instructions into a dictionary for O(1) lookup.
    # Format: {(state, symbol): (new_symbol, direction, new_state)}
    instructions = {}
    for line in instructions_str:
        parts = line.split()
        key = (parts[0], parts[1])
        value = (parts[2], parts[3], parts[4])
        instructions[key] = value

    # Initialize tape, head position, and current state from starting conditions.
    head_pos = start_tape_str.find('H')
    tape = list(start_tape_str.replace('H', ''))
    current_state = start_state

    # List to store the history of tape configurations for the final output.
    history = []

    # 2. Simulation Loop
    while True:
        # Before each step, ensure the head is within the tape list's bounds.
        # Expand the tape with blanks '_' if the head has moved off either end.
        if head_pos < 0:
            tape.insert(0, '_')
            head_pos = 0
        elif head_pos >= len(tape):
            tape.append('_')
        
        # Format the current tape configuration string.
        # e.g., tape=['1', '1', '0'], head_pos=1 -> "1H10"
        current_tape_view = list(tape)
        current_tape_view.insert(head_pos, 'H')
        history.append("".join(current_tape_view))

        # Check for the halting condition.
        if current_state == 'halt':
            break

        # Read the symbol under the head.
        symbol_under_head = tape[head_pos]

        # Find the instruction for the current state and symbol.
        action = instructions[(current_state, symbol_under_head)]
        new_symbol, direction, new_state = action

        # Execute the instruction.
        # a) Write the new symbol onto the tape.
        tape[head_pos] = new_symbol
        
        # b) Move the head left or right.
        if direction == 'r':
            head_pos += 1
        elif direction == 'l':
            head_pos -= 1
        
        # c) Transition to the new state.
        current_state = new_state

    # 3. Output
    # Print the comma-separated history of all tape states.
    print(','.join(history))

simulate_turing_machine()
<<<1H10,10H0,100H_,10H0_,1H00_,H100_,H_100_,_H100_>>>