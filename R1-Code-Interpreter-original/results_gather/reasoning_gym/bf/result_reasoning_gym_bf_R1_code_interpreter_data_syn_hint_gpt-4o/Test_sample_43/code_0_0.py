# Simulating the Brainf*ck code execution
def brainfuck_interpreter():
    tape = [0] * 30000
    pointer = 0
    output = []

    # Initial setup
    tape[pointer] = 0
    pointer += 1
    tape[pointer] = 0
    pointer -= 1

    # Increment the first cell by 8
    tape[pointer] += 8

    # Loop to multiply by 14
    while tape[pointer] != 0:
        tape[pointer] -= 1
        pointer += 1
        tape[pointer] += 14
        pointer -= 1

    # Move to the first cell, output, and decrement
    pointer -= 1
    output.append(chr(tape[pointer]))
    tape[pointer] -= 1

    # Increment by 6 and output
    tape[pointer] += 6
    output.append(chr(tape[pointer]))
    tape[pointer] -= 1

    # Increment by 5 and output
    tape[pointer] += 5
    output.append(chr(tape[pointer]))

    return ''.join(output)

print(brainfuck_interpreter())