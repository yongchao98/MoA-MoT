import sys

def solve_turing_machine():
    """
    Simulates a Turing machine based on a given program and initial tape state.
    """

    # --- Configuration from the problem description ---
    program_str = """
    0 1 0 r 0
    0 0 0 r 0
    0 _ _ l 1
    1 0 0 l 1
    1 1 1 l 1
    1 _ _ r halt
    """
    initial_tape_str = "1H10"
    initial_state_num = 0

    # Step 1: Parse the program into a dictionary for easy lookup
    instructions = {}
    for line in program_str.strip().split('\n'):
        parts = line.strip().split()
        if len(parts) == 5:
            current_st, current_sym, new_sym, direction, new_st_str = parts
            
            # Convert state to int, but keep 'halt' as a string
            current_st_key = int(current_st)
            new_st_val = int(new_st_str) if new_st_str != 'halt' else 'halt'
                
            instructions[(current_st_key, current_sym)] = (new_sym, direction, new_st_val)

    # Step 2: Initialize Tape, Head, and State
    head_start_index = initial_tape_str.find('H')
    tape = list(initial_tape_str.replace('H', ''))
    head_pos = head_start_index
    current_state = initial_state_num
    
    history = []

    # Helper function to format the tape state for output
    def format_tape_state(current_tape, current_head_pos):
        if current_head_pos < 0:
            # Head is on an implicit blank to the left
            return 'H_' + "".join(current_tape)
        elif current_head_pos >= len(current_tape):
            # Head is on an implicit blank to the right
            return "".join(current_tape) + 'H_'
        else:
            # Head is on a symbol within the tape
            prefix = "".join(current_tape[0:current_head_pos])
            symbol_at_head = current_tape[current_head_pos]
            suffix = "".join(current_tape[current_head_pos + 1:])
            return prefix + 'H' + symbol_at_head + suffix

    # Step 3: Run the simulation loop
    while current_state != 'halt':
        # Add current configuration to history
        history.append(format_tape_state(tape, head_pos))

        # Read symbol at the head position
        current_symbol = '_'
        if 0 <= head_pos < len(tape):
            current_symbol = tape[head_pos]

        # Get the corresponding instruction
        instruction = instructions.get((current_state, current_symbol))
        if not instruction:
            # This case assumes a well-defined program where this doesn't happen
            break
        
        new_symbol, direction, new_state = instruction
        
        # Write the new symbol to the tape, extending it if necessary
        if head_pos < 0:
            tape.insert(0, new_symbol)
            head_pos = 0 # Head is now at the new first element
        elif head_pos >= len(tape):
            tape.append(new_symbol)
            head_pos = len(tape) - 1 # Head is now at the new last element
        else:
            tape[head_pos] = new_symbol
            
        # Move the head
        if direction == 'r':
            head_pos += 1
        elif direction == 'l':
            head_pos -= 1
            
        # Transition to the new state
        current_state = new_state

    # Add the final halting configuration to the history
    history.append(format_tape_state(tape, head_pos))

    # Step 4: Print the final result
    print(",".join(history))

solve_turing_machine()