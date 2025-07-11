import collections

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
    initial_state = "0"

    # 1. Parse instructions into a more usable format
    rules = collections.defaultdict(dict)
    for line in instructions_str.strip().split('\n'):
        parts = line.strip().split()
        if len(parts) == 5:
            state, symbol, new_symbol, direction, new_state = parts
            rules[state][symbol] = (new_symbol, direction, new_state)

    # 2. Initialize the tape and head
    head_pos = initial_tape_str.find('H')
    tape_list = list(initial_tape_str.replace('H', ''))
    
    current_state = initial_state
    history = []

    def format_tape_state(tape, head):
        """Formats the tape and head position into the required string format."""
        tape_str = "".join(tape)
        if head < 0:
            # Head is to the left of the tape
            return f"_H{tape_str}"
        elif head >= len(tape):
            # Head is to the right of the tape
            return f"{tape_str}H_"
        else:
            # Head is on the tape
            return f"{tape_str[:head]}H{tape_str[head:]}"

    # 3. Simulation loop
    while current_state != "halt":
        # Add current configuration to history
        history.append(format_tape_state(tape_list, head_pos))

        # Read the symbol at the head's position
        current_symbol = '_'
        if 0 <= head_pos < len(tape_list):
            current_symbol = tape_list[head_pos]

        # Get the instruction
        new_symbol, direction, new_state = rules[current_state][current_symbol]
        
        # Write the new symbol to the tape
        # Note: We only modify the tape for non-blank symbols, or if the head is off-tape
        # and we need to expand it. However, the given rules don't write non-blanks off-tape.
        if new_symbol != '_':
            if head_pos < 0:
                # This case is not hit by the given rules, but included for completeness.
                # Writing 'X' left of tape ['a','b'] makes it ['X','a','b'] with head at new 0.
                tape_list.insert(0, new_symbol)
                head_pos = 0 
            elif head_pos >= len(tape_list):
                 # This case is not hit by the given rules.
                tape_list.append(new_symbol)
            else:
                tape_list[head_pos] = new_symbol
        
        # Move the head
        if direction == 'r':
            head_pos += 1
        elif direction == 'l':
            head_pos -= 1

        # Transition to the new state
        current_state = new_state
        
    # Add the final halting state to history
    history.append(format_tape_state(tape_list, head_pos))

    # 4. Print the final result
    print(",".join(history))

solve()
<<<1H10,10H0,100H_,10H0,1H00,H100,_H100,H100>>>