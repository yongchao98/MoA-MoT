def solve_turing_machine():
    """
    Simulates a Turing machine based on a given program and initial tape state.
    """
    # 1. Parse Inputs
    instructions_str = """0 1 0 r 0
0 0 0 r 0
0 _ _ l 1
1 0 0 l 1
1 1 1 l 1
1 _ _ r halt"""
    initial_tape_str = "1H10"
    initial_state = "0"

    instructions = {}
    for line in instructions_str.strip().split('\n'):
        parts = line.split()
        key = (parts[0], parts[1])
        value = (parts[2], parts[3], parts[4])
        instructions[key] = value

    head_pos = initial_tape_str.find('H')
    tape = list(initial_tape_str.replace('H', ''))

    # 2. Initialize State
    current_state = initial_state
    history = [initial_tape_str]

    # 3. Create a Formatting Function
    def format_tape_state(current_tape, current_head_pos):
        """
        Formats the tape and head position into the required string representation,
        trimming blanks unless the head is on them.
        """
        # Find the first and last non-blank characters to determine content boundaries
        first_non_blank = -1
        last_non_blank = -1
        for i, char in enumerate(current_tape):
            if char != '_':
                if first_non_blank == -1:
                    first_non_blank = i
                last_non_blank = i

        # If the tape is all blanks, the head is on a blank
        if first_non_blank == -1:
            return "H_"

        # The view must include all non-blank content and the head
        start_idx = min(first_non_blank, current_head_pos)
        end_idx = max(last_non_blank, current_head_pos)

        # Build the string representation for the relevant tape slice
        output_parts = []
        for i in range(start_idx, end_idx + 1):
            if i == current_head_pos:
                output_parts.append('H')
            output_parts.append(current_tape[i])

        return "".join(output_parts)

    # 4. Run Simulation Loop
    while current_state != 'halt':
        # Handle Tape Boundaries (extend tape if head moves off the edge)
        if head_pos < 0:
            tape.insert(0, '_')
            head_pos = 0  # Head is now at the new start
        elif head_pos >= len(tape):
            tape.append('_')

        # Get current symbol and find the instruction
        current_symbol = tape[head_pos]
        rule_key = (current_state, current_symbol)
        
        if rule_key not in instructions:
            print(f"Error: No instruction found for state '{current_state}' and symbol '{current_symbol}'")
            break
            
        new_symbol, direction, new_state = instructions[rule_key]

        # Execute instruction
        tape[head_pos] = new_symbol
        if direction == 'r':
            head_pos += 1
        elif direction == 'l':
            head_pos -= 1
        current_state = new_state

        # Record the new state
        history.append(format_tape_state(tape, head_pos))

    # 5. Final Output
    print(",".join(history))

solve_turing_machine()
<<<1H10,10H0,100H_,10H0,1H00,H100,H_100,_H100>>>