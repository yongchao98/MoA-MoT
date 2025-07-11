import collections

def solve():
    """
    Simulates the execution of a Turing machine and prints the history of its states.
    """
    # Define the Turing machine's program and initial state
    instructions_str = """
    0 1 0 r 0
    0 0 0 r 0
    0 _ _ l 1
    1 0 0 l 1
    1 1 1 l 1
    1 _ _ r halt
    """
    initial_tape_str = "1H10"
    initial_state = "0"

    # Parse the instructions into a dictionary for easy lookup
    instructions = {}
    for line in instructions_str.strip().split('\n'):
        parts = line.strip().split()
        # Key: (current_state, current_symbol), Value: (new_symbol, direction, new_state)
        instructions[(parts[0], parts[1])] = (parts[2], parts[3], parts[4])

    # Parse the initial tape string to initialize the tape and head position
    head_pos = initial_tape_str.find('H')
    tape_list = list(initial_tape_str.replace('H', ''))

    # Initialize machine state
    current_state = initial_state
    history = []

    # Main simulation loop
    while current_state != "halt":
        # --- Format and store the current state ---
        # Find the first and last non-blank symbols to correctly represent the tape
        first_non_blank = -1
        last_non_blank = -1
        for i, char in enumerate(tape_list):
            if char != '_':
                if first_non_blank == -1:
                    first_non_blank = i
                last_non_blank = i
        
        formatted_tape = ""
        if first_non_blank == -1: # Tape is all blanks
             formatted_tape = "H_"
        else:
            # Check if head is outside the active area
            if head_pos < first_non_blank:
                active_part = "".join(tape_list[first_non_blank:last_non_blank+1])
                formatted_tape = f"H_{active_part}"
            elif head_pos > last_non_blank:
                active_part = "".join(tape_list[first_non_blank:last_non_blank+1])
                formatted_tape = f"{active_part}H_"
            else: # Head is on the active area
                active_part_list = tape_list[first_non_blank:last_non_blank+1]
                relative_head = head_pos - first_non_blank
                active_part_list.insert(relative_head, 'H')
                formatted_tape = "".join(active_part_list)
        history.append(formatted_tape)
        
        # --- Execute one step of the Turing machine ---
        
        # Read the current symbol
        current_symbol = '_'
        if 0 <= head_pos < len(tape_list):
            current_symbol = tape_list[head_pos]

        # Find the instruction to execute
        rule = instructions.get((current_state, current_symbol))
        if not rule:
            # This case should not happen in a well-defined machine
            break
        new_symbol, direction, new_state = rule

        # Write the new symbol, extending tape if necessary
        if head_pos < 0:
            tape_list.insert(0, new_symbol)
            head_pos = 0
        elif head_pos >= len(tape_list):
            tape_list.append(new_symbol)
        else:
            tape_list[head_pos] = new_symbol

        # Move the head
        if direction == 'r':
            head_pos += 1
        elif direction == 'l':
            head_pos -= 1

        # Update the state
        current_state = new_state

    # Add the final halting state to the history
    # (The formatting logic is the same as inside the loop)
    first_non_blank = -1
    last_non_blank = -1
    for i, char in enumerate(tape_list):
        if char != '_':
            if first_non_blank == -1:
                first_non_blank = i
            last_non_blank = i
    
    formatted_tape = ""
    if first_non_blank == -1:
         formatted_tape = "H_"
    else:
        if head_pos < first_non_blank:
            active_part = "".join(tape_list[first_non_blank:last_non_blank+1])
            formatted_tape = f"H_{active_part}"
        elif head_pos > last_non_blank:
            active_part = "".join(tape_list[first_non_blank:last_non_blank+1])
            formatted_tape = f"{active_part}H_"
        else:
            active_part_list = tape_list[first_non_blank:last_non_blank+1]
            relative_head = head_pos - first_non_blank
            active_part_list.insert(relative_head, 'H')
            formatted_tape = "".join(active_part_list)
    history.append(formatted_tape)

    # Print the final result
    print(",".join(history))

solve()