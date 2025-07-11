def solve_turing_machine():
    """
    Simulates a Turing machine based on a given program and initial tape state.
    """
    # 1. Parse Inputs
    program_str = """
    0 1 0 r 0
    0 0 0 r 0
    0 _ _ l 1
    1 0 0 l 1
    1 1 1 l 1
    1 _ _ r halt
    """
    instructions = {}
    for line in program_str.strip().split('\n'):
        parts = line.strip().split()
        # Key: (current_state, current_symbol), Value: (new_symbol, direction, new_state)
        instructions[(parts[0], parts[1])] = (parts[2], parts[3], parts[4])

    initial_tape_str = "1H10"
    initial_state = '0'

    # 2. Initialize
    head_pos = initial_tape_str.find('H')
    tape = list(initial_tape_str.replace('H', ''))
    current_state = initial_state
    history = [initial_tape_str]

    def format_tape(current_tape, current_head_pos):
        """
        Formats the tape and head position into the required string representation,
        trimming unnecessary blank symbols.
        """
        # Create a string representation with 'H' inserted before the symbol at the head
        s = ""
        if current_head_pos < 0:
            # Head moved left off the tape, implying a blank
            s = "H_" + "".join(current_tape)
        elif current_head_pos >= len(current_tape):
            # Head moved right off the tape, implying a blank
            s = "".join(current_tape) + "H_"
        else:
            s = "".join(current_tape[:current_head_pos]) + "H" + "".join(current_tape[current_head_pos:])

        # Trim leading/trailing blanks, but not if the head is adjacent to them
        start = 0
        end = len(s)
        
        # Trim leading blanks: move 'start' forward if it's a blank
        # and the next character isn't the head.
        while start < len(s) - 1 and s[start] == '_' and s[start+1] != 'H':
            start += 1
            
        # Trim trailing blanks: move 'end' backward if it's a blank
        # and the previous character isn't the head.
        while end > 1 and s[end-1] == '_' and s[end-2] != 'H':
            end -= 1
            
        return s[start:end]

    # 3. Simulation Loop
    while current_state != 'halt':
        # If head is off the tape, expand the tape with a blank symbol
        if head_pos < 0:
            tape.insert(0, '_')
            head_pos = 0  # Head is now at the new start
        elif head_pos >= len(tape):
            tape.append('_')

        # Read symbol and find the corresponding instruction
        current_symbol = tape[head_pos]
        key = (current_state, current_symbol)
        
        if key not in instructions:
            print(f"Error: No instruction found for state '{current_state}' and symbol '{current_symbol}'")
            break
            
        new_symbol, direction, new_state = instructions[key]

        # Execute the instruction
        tape[head_pos] = new_symbol
        
        if direction == 'r':
            head_pos += 1
        elif direction == 'l':
            head_pos -= 1
            
        current_state = new_state

        # Format the new state and add to history
        formatted_state = format_tape(tape, head_pos)
        history.append(formatted_state)

    # 4. Final Output
    print(",".join(history))

solve_turing_machine()
<<<1H10,10H0,100H_,10H0,1H00,H100,H_100,_H100>>>