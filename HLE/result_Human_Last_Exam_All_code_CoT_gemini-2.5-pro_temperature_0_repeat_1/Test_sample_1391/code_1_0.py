import collections

def simulate_turing_machine():
    """
    Simulates a Turing machine based on a predefined program and initial tape state.
    """
    # Define the Turing machine's program
    instructions_str = """
    0 1 0 r 0
    0 0 0 r 0
    0 _ _ l 1
    1 0 0 l 1
    1 1 1 l 1
    1 _ _ r halt
    """
    
    # Define the initial state
    initial_tape_str = "1H10"
    initial_state = '0'

    # Parse the program instructions into a dictionary
    program = {}
    for line in instructions_str.strip().split('\n'):
        parts = line.strip().split()
        program[(parts[0], parts[1])] = (parts[2], parts[3], parts[4])

    # Initialize the tape and head position
    head_pos = initial_tape_str.find('H')
    tape_list = list(initial_tape_str.replace('H', ''))
    
    # Use a deque for efficient appends and prepends
    tape = collections.deque(tape_list)
    
    current_state = initial_state
    history = []

    def format_tape_state(tape_deque, head_index):
        """Formats the tape and head position into a string."""
        temp_tape = list(tape_deque)
        
        # Handle head being off the left end
        if head_index < 0:
            # The head is at some position to the left of the tape
            # The string should be H_...
            return 'H_' + "".join(temp_tape)
            
        # Handle head being off the right end
        if head_index >= len(temp_tape):
            # The head is at the position right after the tape ends
            # The string should be ...H_
            return "".join(temp_tape) + 'H_'
            
        # Head is on the tape
        temp_tape.insert(head_index, 'H')
        return "".join(temp_tape)

    def format_tape_state_for_history(tape_deque, head_index):
        """
        A more robust formatter that shows the tape content around the head,
        including blanks if the head is on them.
        """
        tape_str = ""
        
        # Create a temporary representation for display
        display_tape = list(tape_deque)
        display_head = head_index
        
        # Adjust for head off the left
        if head_index < 0:
            padding = ['_'] * (-head_index)
            display_tape = padding + display_tape
            display_head = 0
        
        # Adjust for head off the right
        if head_index >= len(display_tape):
            padding = ['_'] * (head_index - len(display_tape) + 1)
            display_tape = display_tape + padding

        # Insert the head marker
        display_tape.insert(display_head, 'H')
        return "".join(display_tape)


    while current_state != 'halt':
        # Record the current state
        # To match the example output format, we need to show the tape's content
        # with blanks if the head is on them.
        
        # Let's trace the tape content to generate the output string
        # This is a bit tricky, so let's build the string manually
        
        temp_tape_list = list(tape)
        output_parts = []
        
        # Part before the head
        if head_pos < 0:
            output_parts.append("H_")
            output_parts.extend(temp_tape_list)
        elif head_pos >= len(temp_tape_list):
            output_parts.extend(temp_tape_list)
            output_parts.append("H_")
        else:
            output_parts.extend(temp_tape_list[0:head_pos])
            output_parts.append("H")
            output_parts.extend(temp_tape_list[head_pos:])
        
        history.append("".join(output_parts))

        # Read the symbol under the head
        if 0 <= head_pos < len(tape):
            current_symbol = tape[head_pos]
        else:
            current_symbol = '_'

        # Get the instruction
        new_symbol, direction, new_state = program[(current_state, current_symbol)]

        # Execute the instruction: Write
        if head_pos < 0:
            # Pad left
            for _ in range(-head_pos):
                tape.appendleft('_')
            head_pos = 0
        elif head_pos >= len(tape):
            # Pad right
            for _ in range(head_pos - len(tape) + 1):
                tape.append('_')
        
        tape[head_pos] = new_symbol

        # Execute the instruction: Move
        if direction == 'r':
            head_pos += 1
        elif direction == 'l':
            head_pos -= 1

        # Update state
        current_state = new_state

    # Add the final halting state to history
    temp_tape_list = list(tape)
    output_parts = []
    if head_pos < 0:
        output_parts.append("H_")
        output_parts.extend(temp_tape_list)
    elif head_pos >= len(temp_tape_list):
        output_parts.extend(temp_tape_list)
        output_parts.append("H_")
    else:
        output_parts.extend(temp_tape_list[0:head_pos])
        output_parts.append("H")
        output_parts.extend(temp_tape_list[head_pos:])
    history.append("".join(output_parts))

    # Print the final result
    print(",".join(history))

simulate_turing_machine()