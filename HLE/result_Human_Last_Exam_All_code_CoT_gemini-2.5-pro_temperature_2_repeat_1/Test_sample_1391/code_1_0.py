import sys
from io import StringIO

def solve():
    """
    Simulates a Turing machine based on a given program and initial tape state.
    """
    
    # 1. Define the program and initial conditions
    instructions_str = """
    0 1 0 r 0
    0 0 0 r 0
    0 _ _ l 1
    1 0 0 l 1
    1 1 1 l 1
    1 _ _ r halt
    """
    initial_tape_with_head = "1H10"
    initial_state = '0'

    # Parse instructions into a dictionary
    program = {}
    for line in instructions_str.strip().split('\n'):
        parts = line.strip().split()
        if len(parts) == 5:
            current_s, current_sym, new_sym, direction, new_s = parts
            program[(current_s, current_sym)] = (new_sym, direction, new_s)

    # Parse initial tape and head position
    head_pos = initial_tape_with_head.find('H')
    tape_list = list(initial_tape_with_head.replace('H', ''))
    
    current_state = initial_state
    head = head_pos
    tape = tape_list
    
    history = []

    def format_state(current_tape, current_head):
        """
        Formats the tape and head position into the required string format.
        It trims leading/trailing blanks unless the head is on them.
        """
        
        # Create a working copy
        tape_copy = list(current_tape)

        # Find the start and end of non-blank symbols
        start_index = 0
        while start_index < len(tape_copy) and tape_copy[start_index] == '_':
            start_index += 1

        end_index = len(tape_copy)
        while end_index > start_index and tape_copy[end_index - 1] == '_':
            end_index -= 1

        # If the tape is effectively all blanks
        if start_index == end_index:
            return "H_"

        # Slice the relevant part of the tape
        sub_tape = tape_copy[start_index:end_index]

        if current_head < start_index:
            # Head is on a leading blank
            return "H_" + "".join(sub_tape)
        elif current_head >= end_index:
            # Head is on a trailing blank
            return "".join(sub_tape) + "H_"
        else:
            # Head is within the non-blank section
            relative_head = current_head - start_index
            return "".join(sub_tape[:relative_head]) + "H" + "".join(sub_tape[relative_head:])

    # Record the initial state
    history.append(initial_tape_with_head)

    # 2. Main simulation loop
    while current_state != 'halt':
        # Determine the symbol under the head
        if 0 <= head < len(tape):
            current_symbol = tape[head]
        else:
            current_symbol = '_'

        # Find the corresponding instruction
        if (current_state, current_symbol) not in program:
            # No instruction found, implies halting
            break
            
        new_symbol, direction, new_state = program[(current_state, current_symbol)]

        # Execute the instruction: write to the tape
        if head < 0:
            # Prepend to tape and adjust head
            tape.insert(0, new_symbol)
            head = 0
        elif head >= len(tape):
            # Pad with blanks if needed and append
            tape.extend(['_'] * (head - len(tape)))
            tape.append(new_symbol)
        else:
            tape[head] = new_symbol

        # Update state and move head
        current_state = new_state
        
        if direction == 'l':
            head -= 1
        elif direction == 'r':
            head += 1
            
        # Record the new state after the move
        history.append(format_state(tape, head))

    # Print the final comma-separated result
    print(",".join(history))

solve()