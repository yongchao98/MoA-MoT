import sys

def solve():
    """
    Simulates a Turing machine based on a given program and initial tape state.
    """
    instructions_str = """
0 1 0 r 0
0 0 0 r 0
0 _ _ l 1
1 0 0 l 1
1 1 1 l 1
1 _ _ r halt
    """
    initial_tape_str = "1H10"
    initial_state = '0'

    # Parse instructions into a dictionary for easy lookup
    program = {}
    for line in instructions_str.strip().split('\n'):
        parts = line.split()
        current_state, current_symbol, new_symbol, direction, new_state = parts
        program[(current_state, current_symbol)] = (new_symbol, direction, new_state)

    # Parse initial tape and head position
    head_pos = initial_tape_str.find('H')
    tape = list(initial_tape_str.replace('H', ''))
    
    state = initial_state
    history = []

    def format_tape_state(tape_list, head):
        """
        Formats the tape and head position into the specified string format.
        It trims leading/trailing blanks unless the head is on them.
        """
        temp_tape = list(tape_list)
        
        # A working copy of the head position for formatting
        format_head = head

        # Temporarily expand tape for formatting if head is out of bounds
        while format_head < 0:
            temp_tape.insert(0, '_')
            format_head += 1
        while format_head >= len(temp_tape):
            temp_tape.append('_')

        # Insert 'H' marker
        res_list = [c for c in temp_tape]
        res_list[format_head] = 'H' + res_list[format_head]

        # Find bounds for trimming
        start_idx, end_idx = 0, len(res_list) - 1

        # Trim leading blanks, stopping if 'H' is encountered
        while start_idx < end_idx and res_list[start_idx] == '_' and 'H' not in res_list[start_idx]:
            start_idx += 1
        
        # Trim trailing blanks, stopping if 'H' is encountered
        while end_idx > start_idx and res_list[end_idx] == '_' and 'H' not in res_list[end_idx]:
            end_idx -= 1
        
        return "".join(res_list[start_idx:end_idx + 1])


    while state != 'halt':
        # 1. Record the current state configuration
        history.append(format_tape_state(tape, head_pos))

        # 2. Read the current symbol
        if 0 <= head_pos < len(tape):
            current_symbol = tape[head_pos]
        else:
            current_symbol = '_'

        # 3. Get the instruction to execute
        key = (state, current_symbol)
        if key not in program:
            # Implicit halt if no instruction found
            break
            
        new_symbol, direction, new_state = program[key]

        # 4. Expand tape if head is outside the current bounds
        while head_pos < 0:
            tape.insert(0, '_')
            head_pos += 1
        while head_pos >= len(tape):
            tape.append('_')

        # 5. Write the new symbol
        tape[head_pos] = new_symbol

        # 6. Move the head
        if direction == 'r':
            head_pos += 1
        elif direction == 'l':
            head_pos -= 1

        # 7. Transition to the new state
        state = new_state
    
    # Record the final halting state
    history.append(format_tape_state(tape, head_pos))

    print(','.join(history))

solve()
<<<1H10,10H0,100H_,10H0,1H00,H100,H_100,_H100>>>