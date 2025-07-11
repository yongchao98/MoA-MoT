import sys

def solve_turing_machine():
    """
    Simulates a Turing machine based on a hardcoded program and initial tape state.
    """
    # 1. Define the Turing Machine Program and Initial State
    instructions = [
        "0 1 0 r 0",
        "0 0 0 r 0",
        "0 _ _ l 1",
        "1 0 0 l 1",
        "1 1 1 l 1",
        "1 _ _ r halt"
    ]
    initial_tape_str = "1H10"
    machine_state = '0'

    # 2. Parse rules into a lookup dictionary
    rules = {}
    for instruction in instructions:
        parts = instruction.split()
        current_state, current_symbol, new_symbol, direction, new_state = parts
        rules[(current_state, current_symbol)] = (new_symbol, direction, new_state)

    # 3. Initialize tape and head position
    head_pos = initial_tape_str.find('H')
    tape = list(initial_tape_str.replace('H', ''))

    history = []

    def format_tape_state(current_tape, current_head_pos):
        """
        Formats the tape and head position into the required string representation.
        It shows the minimal part of the tape containing all non-blank symbols and the head.
        """
        non_blank_indices = [i for i, char in enumerate(current_tape) if char != '_']
        
        # Determine the slice of the tape to display
        if not non_blank_indices:
            # If tape is all blanks, just show the head on a single blank
            min_idx = current_head_pos
            max_idx = current_head_pos
        else:
            indices_to_consider = non_blank_indices + [current_head_pos]
            min_idx = min(indices_to_consider)
            max_idx = max(indices_to_consider)
        
        display_tape = current_tape[min_idx : max_idx + 1]
        display_head_pos = current_head_pos - min_idx
        
        # Insert 'H' to mark the head position
        final_tape_list = list(display_tape)
        final_tape_list.insert(display_head_pos, 'H')
        
        return "".join(final_tape_list)

    # 4. Simulation loop
    while machine_state != 'halt':
        # Record current state
        history.append(format_tape_state(tape, head_pos))

        # Read the symbol under the head
        current_symbol = tape[head_pos]

        # Find and apply the rule
        rule = rules.get((machine_state, current_symbol))
        if not rule:
            # If no rule matches, halt implicitly
            break

        new_symbol, direction, new_state = rule

        # Write the new symbol to the tape
        tape[head_pos] = new_symbol
        
        # Move the head
        if direction == 'r':
            head_pos += 1
        elif direction == 'l':
            head_pos -= 1

        # Update the machine state
        machine_state = new_state
        
        # Expand tape if head moves off the ends
        if head_pos >= len(tape):
            tape.append('_')
        elif head_pos < 0:
            tape.insert(0, '_')
            head_pos = 0  # After prepending, new head position is 0

    # 5. Record the final halting state
    history.append(format_tape_state(tape, head_pos))

    # 6. Print the comma-separated history
    print(",".join(history))

solve_turing_machine()