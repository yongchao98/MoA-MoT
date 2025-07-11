import sys

def solve():
    """
    Simulates a Turing machine and prints the history of its tape states.
    """
    program_text = """
0 1 0 r 0
0 0 0 r 0
0 _ _ l 1
1 0 0 l 1
1 1 1 l 1
1 _ _ r halt
"""
    initial_tape_str = "1H10"
    initial_state = '0'

    # Parse the program instructions into a dictionary
    program = {}
    for line in program_text.strip().split('\n'):
        parts = line.split()
        if len(parts) == 5:
            current_state, current_symbol, new_symbol, direction, new_state = parts
            program[(current_state, current_symbol)] = (new_symbol, direction, new_state)

    # Parse the initial tape configuration
    h_pos = initial_tape_str.find('H')
    tape = list(initial_tape_str.replace('H', ''))
    head = h_pos
    state = initial_state
    
    history = []

    def format_tape_state(current_tape, current_head):
        """
        Formats the tape list and head position into the required string representation,
        trimming unnecessary blanks.
        """
        display_tape = list(current_tape)
        display_head = current_head

        # Trim leading blanks, but only if the head is not on them
        while len(display_tape) > 0 and display_tape[0] == '_' and display_head > 0:
            display_tape.pop(0)
            display_head -= 1

        # Trim trailing blanks, but only if the head is not on them
        while len(display_tape) > 0 and display_tape[-1] == '_' and display_head < len(display_tape) - 1:
            display_tape.pop()
        
        # If the tape is now empty (all blanks), the representation is H_
        if not display_tape:
            return "H_"
            
        # Insert the head 'H' into the string
        if 0 <= display_head < len(display_tape):
             display_tape.insert(display_head, 'H')
        elif display_head >= len(display_tape): # Head is on a blank to the right
             display_tape.append('H')
             display_tape.append('_')
        # Note: A negative head position is handled by the normalization step before this function is called.

        return "".join(display_tape)

    # Record the initial state
    history.append(format_tape_state(tape, head))

    # Simulation loop
    while state != 'halt':
        # Read the symbol under the head, handling off-tape positions
        current_symbol = '_'
        if 0 <= head < len(tape):
            current_symbol = tape[head]

        # Get the instruction
        instruction = program.get((state, current_symbol))
        if instruction is None:
            # No instruction found, implies an unhandled state or halt
            break
        
        new_symbol, direction, new_state = instruction

        # Write the new symbol to the tape
        if 0 <= head < len(tape):
            tape[head] = new_symbol
        # Off-tape writes will be handled by the normalization after the move

        # Update state
        state = new_state

        # Move head
        if direction == 'r':
            head += 1
        elif direction == 'l':
            head -= 1
        
        # Normalize tape and head for off-tape moves
        if head < 0:
            tape.insert(0, '_')
            head = 0
        elif head >= len(tape):
            # The symbol to write on the new blank space is the one from the instruction
            tape.append(new_symbol)


        # Record the new state after the entire instruction is executed
        history.append(format_tape_state(tape, head))
    
    # Print the final comma-separated list
    final_output = ",".join(history)
    print(final_output)
    print(f"<<<{final_output}>>>")

solve()