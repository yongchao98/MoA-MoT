def solve_turing_machine():
    """
    Simulates a Turing machine based on a given program and initial tape state.
    """

    # 1. Define the machine's program and initial state
    program_str = """
0 1 0 r 0
0 0 0 r 0
0 _ _ l 1
1 0 0 l 1
1 1 1 l 1
1 _ _ r halt
    """
    initial_tape_str = "1H10"
    current_state = "0"

    # 2. Parse the program instructions into a dictionary for easy lookup
    rules = {}
    for line in program_str.strip().split('\n'):
        parts = line.split()
        # Key: (state, symbol), Value: (new_symbol, direction, new_state)
        rules[(parts[0], parts[1])] = (parts[2], parts[3], parts[4])

    # 3. Parse the initial tape configuration
    head_pos = initial_tape_str.find('H')
    tape = list(initial_tape_str.replace('H', ''))

    # List to store the history of tape configurations
    history = []

    def format_tape_state(tape_list, head_idx):
        """
        Formats the tape list and head position into the required string representation.
        It handles tape expansion and trimming of unnecessary blanks.
        """
        temp_tape = list(tape_list)  # Work on a copy

        # Temporarily expand the tape to place the head if it's out of bounds
        while head_idx < 0:
            temp_tape.insert(0, '_')
            head_idx += 1
        while head_idx >= len(temp_tape):
            temp_tape.append('_')

        # Insert the head marker 'H'
        display_list = list(temp_tape)
        display_list.insert(head_idx, 'H')
        s = "".join(display_list)

        # Trim leading blanks, but not if the head is on the first one
        while s.startswith('_') and s[1] != 'H':
            s = s[1:]
        # Trim trailing blanks, but not if the head is on the last one
        while s.endswith('_') and s[-2] != 'H':
            s = s[:-1]
        return s

    # 4. Main simulation loop
    while current_state != "halt":
        # Add the current formatted state to history
        history.append(format_tape_state(tape, head_pos))

        # Read the current symbol, expanding the tape if the head is out of bounds
        current_symbol = '_'
        if head_pos < 0:
            tape.insert(0, '_')
            head_pos = 0
        elif head_pos >= len(tape):
            tape.append('_')
        
        current_symbol = tape[head_pos]

        # Find the rule for the current state and symbol
        rule = rules[(current_state, current_symbol)]
        new_symbol, direction, new_state = rule

        # Execute the rule: Write, Move, Change State
        tape[head_pos] = new_symbol

        if direction == 'l':
            head_pos -= 1
        elif direction == 'r':
            head_pos += 1

        current_state = new_state

    # 5. Add the final halting state to the history
    history.append(format_tape_state(tape, head_pos))

    # 6. Print the comma-separated history
    print(",".join(history))


solve_turing_machine()