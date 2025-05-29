# Simulating the Brainf*ck program
memory = [0] * 30000
pointer = 0

# Executing the program
memory[pointer + 1] = 0
memory[pointer + 2] = 0
pointer += 1
memory[pointer] += 8

while memory[pointer] != 0:
    pointer -= 1
    memory[pointer] += 14
    pointer += 1
    memory[pointer] -= 1

pointer -= 1
output = chr(memory[pointer])
memory[pointer] -= 11
output += chr(memory[pointer])
memory[pointer] += 9
output += chr(memory[pointer])
output += chr(memory[pointer])
memory[pointer] -= 9
output += chr(memory[pointer])

print(output)