import sys

def simulate_turing_machine():
    """
    Simulates a Turing machine based on a defined program and initial tape state.
    """
    # Define the program instructions as a multi-line string
    instructions_str = """
    0 1 0 r 0
    0 0 0 r 0
    0 _ _ l 1
    1 0 0 l 1
    1 1 1 l 1
    1 _ _ r halt
    """
    
    # Parse the instructions into a dictionary for easy lookup
    program = {}
    for line in instructions_str.strip().split('\n'):
        parts = line.strip().split()
        if len(parts) == 5:
            current_state, current_symbol, new_symbol, direction, new_state = parts
            program[(current_state, current_symbol)] = (new_symbol, direction, new_state)

    # Set up the initial configuration of the machine
    initial_tape_str = "1H10"
    initial_state = '0'
    
    # Parse the initial tape string to set the tape content and head position
    head_pos = initial_tape_str.find('H')
    tape_list = list(initial_tape_str.replace('H', ''))
    
    # A dictionary is used to represent the tape to handle its infinite and sparse nature
    tape = {i: symbol for i, symbol in enumerate(tape_list)}

    # Track the boundaries of the tape that have been accessed for printing
    min_tape_idx = 0
    max_tape_idx = len(tape_list) - 1

    current_state = initial_state
    history = []

    # Main simulation loop, continues until the machine halts
    while current_state != 'halt':
        # Update the boundaries of the tape based on the head's position
        min_tape_idx = min(min_tape_idx, head_pos)
        max_tape_idx = max(max_tape_idx, head_pos)
        
        # Format the current state of the tape for output
        tape_str_builder = []
        for i in range(min_tape_idx, max_tape_idx + 1):
            if i == head_pos:
                tape_str_builder.append("H")
            tape_str_builder.append(tape.get(i, '_'))
        
        # This handles the case where the head has moved past the rightmost boundary
        if head_pos > max_tape_idx:
            tape_str_builder.append("H" + tape.get(head_pos, '_'))

        history.append("".join(tape_str_builder))

        # Determine the next action based on the current state and symbol under the head
        current_symbol = tape.get(head_pos, '_')
        instruction_key = (current_state, current_symbol)
        
        if instruction_key not in program:
            print(f"Error: No instruction found for state='{current_state}', symbol='{current_symbol}'", file=sys.stderr)
            break
            
        new_symbol, direction, new_state = program[instruction_key]

        # Execute the instruction: write to tape, move head, and change state
        tape[head_pos] = new_symbol
        
        if direction == 'r':
            head_pos += 1
        elif direction == 'l':
            head_pos -= 1

        current_state = new_state

    # After halting, record the final state of the tape
    min_tape_idx = min(min_tape_idx, head_pos)
    max_tape_idx = max(max_tape_idx, head_pos)
    
    final_tape_str_builder = []
    for i in range(min_tape_idx, max_tape_idx + 1):
        if i == head_pos:
            final_tape_str_builder.append("H")
        final_tape_str_builder.append(tape.get(i, '_'))

    history.append("".join(final_tape_str_builder))
    
    # Print the complete history as a single comma-separated string
    print(",".join(history))

simulate_turing_machine()
<<<1H10,10H0,100H_,10H0_,1H00_,H100_,H_100_,_H100_>>>