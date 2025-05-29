# Simulating the Brainf*ck program
tape = [0] * 30000
pointer = 0

# Executing the program
tape[pointer + 1] = 10
while tape[pointer + 1] != 0:
    tape[pointer] += 11
    tape[pointer + 1] -= 1

pointer -= 1
output = []
output.append(chr(tape[pointer]))  # 'n'
tape[pointer] += 1
output.append(chr(tape[pointer]))  # 'o'
tape[pointer] += 3
output.append(chr(tape[pointer]))  # 'r'
tape[pointer] += 2
output.append(chr(tape[pointer]))  # 't'
tape[pointer] -= 12
output.append(chr(tape[pointer]))  # 'h'

print(''.join(output))