# Simulating the Brainf*ck program
memory = [0] * 30000
pointer = 0

# Program execution
memory[pointer] = 0
pointer += 1
memory[pointer] = 0
pointer -= 1
pointer += 1
memory[pointer] += 9
while memory[pointer] != 0:
    pointer -= 1
    memory[pointer] += 13
    pointer += 1
    memory[pointer] -= 1
pointer -= 1
output = chr(memory[pointer])
memory[pointer] -= 7
output += chr(memory[pointer])
memory[pointer] += 6
output += chr(memory[pointer])
memory[pointer] -= 11
output += chr(memory[pointer])
memory[pointer] += 3
output += chr(memory[pointer])

print(output)