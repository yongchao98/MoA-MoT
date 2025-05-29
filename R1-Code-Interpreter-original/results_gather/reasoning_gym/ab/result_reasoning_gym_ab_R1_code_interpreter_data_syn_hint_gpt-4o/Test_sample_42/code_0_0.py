def apply_rules(program):
    i = 0
    while i < len(program) - 1:
        if program[i] == "A#" and program[i + 1] == "#A":
            # A# #A becomes nothing
            program.pop(i)
            program.pop(i)
        elif program[i] == "A#" and program[i + 1] == "#B":
            # A# #B becomes #B A#
            program[i], program[i + 1] = program[i + 1], program[i]
            i += 1
        elif program[i] == "B#" and program[i + 1] == "#A":
            # B# #A becomes #A B#
            program[i], program[i + 1] = program[i + 1], program[i]
            i += 1
        elif program[i] == "B#" and program[i + 1] == "#B":
            # B# #B becomes nothing
            program.pop(i)
            program.pop(i)
        else:
            i += 1
    return program

# Initial program
program = ["A#", "#A", "A#", "#A", "B#", "B#", "#A", "A#", "#B", "#B"]

# Apply rules until no more changes can be made
previous_length = -1
while previous_length != len(program):
    previous_length = len(program)
    program = apply_rules(program)

# Print the final state of the program
print(" ".join(program))