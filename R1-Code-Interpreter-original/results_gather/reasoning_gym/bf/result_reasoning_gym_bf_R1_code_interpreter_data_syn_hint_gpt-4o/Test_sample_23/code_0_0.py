# Simulating the Brainf*ck program
def brainfuck_interpreter():
    tape = [0] * 30000
    pointer = 0
    output = []

    # Initial setup
    pointer += 1  # >
    tape[pointer] = 0  # [-]
    pointer += 1  # >
    tape[pointer] = 0  # [-]
    pointer -= 1  # <
    pointer -= 1  # <
    tape[pointer] += 8  # ++++++++
    
    # Loop [<++++++++++++++>-]
    while tape[pointer] != 0:
        pointer -= 1  # <
        tape[pointer] += 14  # ++++++++++++++
        pointer += 1  # >
        tape[pointer] -= 1  # -
    
    pointer -= 1  # <
    output.append(chr(tape[pointer]))  # .
    tape[pointer] += 2  # ++
    output.append(chr(tape[pointer]))  # .
    tape[pointer] -= 9  # ---------
    output.append(chr(tape[pointer]))  # .
    tape[pointer] += 4  # ++++
    output.append(chr(tape[pointer]))  # .
    tape[pointer] += 2  # ++
    output.append(chr(tape[pointer]))  # .
    pointer -= 1  # <

    return ''.join(output)

print(brainfuck_interpreter())