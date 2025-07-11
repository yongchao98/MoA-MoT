import sys

def solve():
    """
    Simulates a Turing machine and prints the history of its states.
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

    # 1. Parse instructions into a dictionary
    rules = {}
    for line in instructions_str.strip().split('\n'):
        parts = line.strip().split()
        current_state, current_symbol, new_symbol, direction, new_state = parts
        rules[(current_state, current_symbol)] = (new_symbol, direction, new_state)

    # 2. Initialize the Turing machine
    head_pos = initial_tape_str.find('H')
    tape_list = list(initial_tape_str.replace('H', ''))
    current_state = initial_state
    history = []

    def format_tape(tape, head):
        """Formats the tape list and head position into a string."""
        return "".join(tape[:head]) + "H" + "".join(tape[head:])

    # 3. Simulation Loop
    while current_state != 'halt':
        # Add current configuration to history
        history.append(format_tape(tape_list, head_pos))

        # Read symbol at head position
        symbol_at_head = tape_list[head_pos]

        # Find and execute the rule
        if (current_state, symbol_at_head) in rules:
            new_symbol, direction, new_state = rules[(current_state, symbol_at_head)]

            # Write the new symbol
            tape_list[head_pos] = new_symbol

            # Move the head
            if direction == 'r':
                head_pos += 1
            elif direction == 'l':
                head_pos -= 1

            # Update the state
            current_state = new_state

            # Extend tape if head is out of bounds
            if head_pos < 0:
                tape_list.insert(0, '_')
                head_pos = 0
            elif head_pos >= len(tape_list):
                tape_list.append('_')
        else:
            # Should not happen with a well-defined machine, but as a safeguard
            break

    # Add the final halting state to history
    history.append(format_tape(tape_list, head_pos))

    # 4. Print the final result
    print(",".join(history))

solve()