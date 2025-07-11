import sys

def format_tape(tape, head):
    """Formats the tape and head position into the required string representation."""
    
    # Handle cases where the head is outside the current tape list
    if head < 0:
        # Create representation like H_...
        # Here we assume head is at -1 relative to the tape start
        full_str = "H_" + "".join(tape)
    elif head >= len(tape):
        # Create representation like ...H_
        # Here we assume head is at len(tape) relative to the tape start
        full_str = "".join(tape) + "H_"
    else:
        # Head is inside the tape
        pre = "".join(tape[:head])
        symbol = tape[head]
        post = "".join(tape[head+1:])
        full_str = f"{pre}H{symbol}{post}"
        
    # Trim leading blanks unless the head is on the first char and it's a blank.
    # The representation for head on a leading blank is like _H... or H_...
    # The check "!full_str.startswith('H')" handles "_H" cases correctly.
    while full_str.startswith('_') and not full_str.startswith('_H'):
        full_str = full_str[1:]
        
    # Trim trailing blanks unless the head is on the last char and it's a blank.
    while full_str.endswith('_') and not full_str.endswith('H_'):
        full_str = full_str[:-1]
        
    return full_str

def solve():
    """
    Simulates the Turing machine as described in the problem.
    """
    instructions_str = [
        "0 1 0 r 0",
        "0 0 0 r 0",
        "0 _ _ l 1",
        "1 0 0 l 1",
        "1 1 1 l 1",
        "1 _ _ r halt"
    ]
    
    initial_tape_str = "1H10"
    initial_state = '0'
    
    # Parse instructions into a dictionary for easy lookup
    rules = {}
    for line in instructions_str:
        parts = line.split()
        rules[(parts[0], parts[1])] = (parts[2], parts[3], parts[4])
        
    # Initialize tape and head
    h_pos = initial_tape_str.find('H')
    tape = list(initial_tape_str.replace('H', ''))
    head = h_pos
    
    state = initial_state
    history = []
    
    # Simulation loop
    while state != 'halt':
        # 1. Record current configuration
        history.append(format_tape(tape, head))
        
        # 2. Read current symbol
        current_symbol = ''
        if head < 0 or head >= len(tape):
            current_symbol = '_'
        else:
            current_symbol = tape[head]
            
        # 3. Get instruction
        new_symbol, direction, new_state = rules[(state, current_symbol)]
        
        # 4. Execute instruction
        # a. Write new symbol
        if head < 0:
            tape.insert(0, new_symbol)
            head = 0  # Head is now at the beginning of the modified tape
        elif head >= len(tape):
            tape.append(new_symbol)
            # Head position remains the same, pointing to the new symbol
        else:
            tape[head] = new_symbol
            
        # b. Move head
        if direction == 'l':
            head -= 1
        elif direction == 'r':
            head += 1
            
        # c. Update state
        state = new_state
        
    # Add the final halting configuration to history
    history.append(format_tape(tape, head))
    
    # Print the comma-separated result
    print(",".join(history))

solve()