def solve_turing_machine():
    """
    Simulates a Turing machine based on a hardcoded program and initial tape state.
    """
    # 1. Define the Turing machine's program and initial configuration.
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

    # 2. Parse the program into a more usable dictionary format.
    # The key is (state, symbol), and the value is (new_symbol, direction, new_state).
    instructions = {}
    for line in program_str.strip().split('\n'):
        parts = line.strip().split()
        if len(parts) == 5:
            instructions[(parts[0], parts[1])] = (parts[2], parts[3], parts[4])

    # 3. Parse the initial tape string to get the tape content and head position.
    h_pos = initial_tape_str.find('H')
    head_pos = h_pos
    tape = list(initial_tape_str.replace('H', ''))

    # 4. Initialize simulation variables.
    current_state = initial_state
    history = []

    # 5. Run the simulation loop until the machine halts.
    while current_state != 'halt':
        # 5a. Format the current tape and head position into a string for the history.
        tape_repr = ""
        if head_pos < 0:
            tape_repr = "H_" + "".join(tape)
        elif head_pos >= len(tape):
            tape_repr = "".join(tape) + "H_"
        else:
            tape_repr = "".join(tape[:head_pos]) + "H" + "".join(tape[head_pos:])
        history.append(tape_repr)

        # 5b. Read the symbol under the head. Handle out-of-bounds as a blank '_'.
        current_symbol = ''
        if 0 <= head_pos < len(tape):
            current_symbol = tape[head_pos]
        else:
            current_symbol = '_'

        # 5c. Find the corresponding instruction.
        instruction = instructions[(current_state, current_symbol)]
        new_symbol, direction, new_state = instruction

        # 5d. Write the new symbol to the tape.
        # If the head is off the tape, expand the tape first.
        if head_pos < 0:
            tape.insert(0, new_symbol)
            head_pos = 0  # Adjust head to point to the new cell at index 0
        elif head_pos >= len(tape):
            tape.append(new_symbol)
        else:
            tape[head_pos] = new_symbol

        # 5e. Update the machine's state.
        current_state = new_state

        # 5f. Move the head.
        if direction == 'r':
            head_pos += 1
        elif direction == 'l':
            head_pos -= 1

    # 6. After halting, record the final state of the tape.
    final_tape_repr = ""
    if head_pos < 0:
        final_tape_repr = "H_" + "".join(tape)
    elif head_pos >= len(tape):
        final_tape_repr = "".join(tape) + "H_"
    else:
        final_tape_repr = "".join(tape[:head_pos]) + "H" + "".join(tape[head_pos:])
    history.append(final_tape_repr)

    # 7. Print the complete history as a comma-separated string.
    print(",".join(history))

solve_turing_machine()
<<<1H10,10H0,100H_,10H0_,1H00_,H100_,H_100_,_H100_>>>