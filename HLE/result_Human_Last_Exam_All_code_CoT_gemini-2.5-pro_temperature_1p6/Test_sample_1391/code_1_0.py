import sys
from io import StringIO

def solve():
    """
    Simulates a Turing machine based on a given program and initial tape state.
    """
    # 1. Define inputs
    instructions_text = """0 1 0 r 0
0 0 0 r 0
0 _ _ l 1
1 0 0 l 1
1 1 1 l 1
1 _ _ r halt
"""
    initial_tape_str = "1H10"
    state = "0"
    blank_symbol = '_'

    # 2. Parse the instructions into a more usable format
    program = {}
    for line in instructions_text.strip().split('\n'):
        parts = line.split()
        if len(parts) == 5:
            current_state, current_symbol, new_symbol, direction, new_state = parts
            program[(current_state, current_symbol)] = (new_symbol, direction, new_state)

    # 3. Parse the initial tape state
    if 'H' in initial_tape_str:
        head_pos = initial_tape_str.find('H')
        tape = list(initial_tape_str.replace('H', ''))
    else:
        head_pos = 0
        tape = list(initial_tape_str)

    # History to store each configuration
    history = []

    def format_tape_state(current_tape, current_head):
        """
        Formats the tape and head position into the required string representation.
        It trims leading/trailing blanks unless the head is positioned on them.
        """
        # Find the indices of the first and last non-blank characters
        first_non_blank = -1
        last_non_blank = -1
        for i, symbol in enumerate(current_tape):
            if symbol != blank_symbol:
                if first_non_blank == -1:
                    first_non_blank = i
                last_non_blank = i
        
        # If the tape is all blanks
        if first_non_blank == -1:
            return f"H{blank_symbol}"
        
        # Determine the slice of the tape to display
        start_index = min(first_non_blank, current_head)
        end_index = max(last_non_blank, current_head)
        
        display_tape = current_tape[start_index : end_index + 1]
        
        # Calculate the head's position relative to the slice
        relative_head_pos = current_head - start_index
        
        # Insert 'H' to mark the head's position
        result = "".join(display_tape[:relative_head_pos]) + 'H' + "".join(display_tape[relative_head_pos:])
        return result

    # 4. Run the simulation loop
    while state != 'halt':
        # Store the current state representation
        history.append(format_tape_state(tape, head_pos))

        # Ensure tape is large enough for the head
        if head_pos < 0:
            tape.insert(0, blank_symbol)
            head_pos = 0  # Head is now at the new start of the tape
        elif head_pos >= len(tape):
            tape.append(blank_symbol)
            
        current_symbol = tape[head_pos]

        # Find the transition rule
        key = (state, current_symbol)
        if key not in program:
            # Implicit halt if no rule is found
            break

        new_symbol, direction, new_state = program[key]

        # Execute the transition
        tape[head_pos] = new_symbol
        if direction == 'l':
            head_pos -= 1
        elif direction == 'r':
            head_pos += 1
        
        state = new_state

    # 5. Add the final halting configuration to the history
    history.append(format_tape_state(tape, head_pos))
    
    # 6. Print the final result
    final_output = ",".join(history)
    print(final_output)
    print(f"<<<{final_output}>>>")

solve()