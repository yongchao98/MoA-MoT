def solve_turing_machine():
    """
    This function simulates a Turing machine based on a hardcoded program and
    initial tape state, printing the history of tape configurations until it halts.
    """
    
    # --- 1. Define the Turing Machine's Program and Initial Conditions ---
    
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

    # --- 2. Initialization ---
    
    # Parse the program string into a dictionary for efficient lookup.
    # The key is a tuple (current_state, current_symbol), and the value is
    # a tuple (new_symbol, direction, new_state).
    instructions = {}
    for line in program_str.strip().split('\n'):
        parts = line.split()
        if len(parts) == 5:
            key = (parts[0], parts[1])
            value = (parts[2], parts[3], parts[4])
            instructions[key] = value

    # Parse the initial tape string to set up the tape and head position.
    # A dictionary represents the infinite tape, mapping integer positions to symbols.
    tape = {}
    tape_symbols_str = initial_tape_str.replace('H', '')
    head_pos = initial_tape_str.find('H')
    
    for i, char in enumerate(tape_symbols_str):
        # We only need to store non-blank symbols.
        if char != '_':
            tape[i] = char

    current_state = initial_state
    history = []

    # --- 3. Helper Function for Formatting ---
    
    def format_tape_state(current_tape, current_head_pos):
        """
        Formats the tape and head position into the required string format.
        It determines the minimal string that includes all non-blank symbols
        and the current head position.
        """
        non_blank_indices = [k for k, v in current_tape.items() if v != '_']
        
        # Determine the range of the tape to display.
        if not non_blank_indices:
            # If the tape is empty, the view is centered on the head.
            min_idx, max_idx = current_head_pos, current_head_pos
        else:
            min_idx = min(min(non_blank_indices), current_head_pos)
            max_idx = max(max(non_blank_indices), current_head_pos)
            
        # Create a list of symbols for the display range.
        display_tape_list = [current_tape.get(i, '_') for i in range(min_idx, max_idx + 1)]
        
        # Calculate the head's position relative to the start of our display string.
        display_head_idx = current_head_pos - min_idx
        
        # Insert 'H' at the correct position.
        display_tape_list.insert(display_head_idx, 'H')
        
        return "".join(display_tape_list)
        
    # --- 4. Simulation Loop ---
    
    # Run the simulation until the machine enters the 'halt' state.
    # A safety break is included to prevent any unforeseen infinite loops.
    for _ in range(100):
        # Record the current state of the tape before executing the instruction.
        history.append(format_tape_state(tape, head_pos))

        # If we have reached the halting state, the simulation is over.
        if current_state == 'halt':
            break

        # Read the symbol under the head. A position not in the tape dict is a blank.
        current_symbol = tape.get(head_pos, '_')
        
        # Find the instruction for the current state and symbol.
        instruction_key = (current_state, current_symbol)
        if instruction_key not in instructions:
            # If no rule is defined, the machine implicitly halts.
            break
            
        new_symbol, direction, new_state = instructions[instruction_key]
        
        # Execute the instruction: write the new symbol to the tape.
        if new_symbol == '_':
            # To keep the tape dictionary sparse, remove keys for blank symbols.
            if head_pos in tape:
                del tape[head_pos]
        else:
            tape[head_pos] = new_symbol
            
        # Move the head left or right.
        if direction == 'r':
            head_pos += 1
        elif direction == 'l':
            head_pos -= 1
            
        # Transition to the new state.
        current_state = new_state
    
    # --- 5. Final Output ---
    # Print the comma-separated history of all tape states.
    print(','.join(history))

solve_turing_machine()
<<<1H10,10H0,100H_,10H0,1H00,H100,H_100,H100>>>