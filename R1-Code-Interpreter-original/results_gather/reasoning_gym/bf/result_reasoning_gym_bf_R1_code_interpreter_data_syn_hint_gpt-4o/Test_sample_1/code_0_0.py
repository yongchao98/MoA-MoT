# Simulating the Brainf*ck program
memory = [0] * 30000
pointer = 0
output = []

# Program execution
memory[pointer] = 0
pointer += 1
memory[pointer] = 0
pointer -= 1
pointer += 1
memory[pointer] += 10
while memory[pointer] != 0:
    pointer -= 1
    memory[pointer] += 10
    pointer += 1
    memory[pointer] -= 1
pointer -= 1
memory[pointer] += 1
output.append(chr(memory[pointer]))
memory[pointer] -= 1
memory[pointer] += 5
output.append(chr(memory[pointer]))
memory[pointer] -= 3
output.append(chr(memory[pointer]))
pointer += 1
memory[pointer] += 3
while memory[pointer] != 0:
    pointer -= 1
    memory[pointer] += 6
    pointer += 1
    memory[pointer] -= 1
pointer -= 1
memory[pointer] += 1
output.append(chr(memory[pointer]))

print(''.join(output))