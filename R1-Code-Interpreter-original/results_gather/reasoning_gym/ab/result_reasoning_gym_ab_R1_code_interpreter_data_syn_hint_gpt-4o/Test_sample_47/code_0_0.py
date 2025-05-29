def apply_rules(program):
    i = 0
    while i < len(program) - 1:
        pair = program[i] + " " + program[i + 1]
        if pair == "A# #A":
            program = program[:i] + program[i+2:]
        elif pair == "A# #B":
            program = program[:i] + ["#B", "A#"] + program[i+2:]
            i += 1
        elif pair == "B# #A":
            program = program[:i] + ["#A", "B#"] + program[i+2:]
            i += 1
        elif pair == "B# #B":
            program = program[:i] + program[i+2:]
        else:
            i += 1
    return program

# Initial program
program = ["B#", "#B", "#A", "B#", "B#", "B#", "#B", "B#", "A#", "#A"]

# Apply rules until no more changes
final_program = apply_rules(program)

# Print the final state of the program
print(final_program)