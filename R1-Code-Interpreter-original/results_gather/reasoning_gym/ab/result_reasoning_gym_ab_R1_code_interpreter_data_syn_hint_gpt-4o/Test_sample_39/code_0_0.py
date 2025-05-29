def compute_program(program):
    i = 0
    while i < len(program) - 1:
        if program[i:i+2] == ['A#', '#A']:
            program = program[:i] + program[i+2:]
            i = 0  # Restart from the beginning
        elif program[i:i+2] == ['A#', '#B']:
            program = program[:i] + ['#B', 'A#'] + program[i+2:]
            i = 0  # Restart from the beginning
        elif program[i:i+2] == ['B#', '#A']:
            program = program[:i] + ['#A', 'B#'] + program[i+2:]
            i = 0  # Restart from the beginning
        elif program[i:i+2] == ['B#', '#B']:
            program = program[:i] + program[i+2:]
            i = 0  # Restart from the beginning
        else:
            i += 1
    return program

initial_program = ['#A', 'A#', '#A', 'B#', 'B#', '#A', 'A#', 'B#', 'A#', 'B#']
final_state = compute_program(initial_program)
print(final_state)