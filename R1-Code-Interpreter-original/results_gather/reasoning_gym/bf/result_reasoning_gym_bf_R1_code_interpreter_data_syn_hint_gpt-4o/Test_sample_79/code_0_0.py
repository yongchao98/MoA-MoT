# Initialize the memory tape and pointer
tape = [0] * 30000
pointer = 0

# Execute the Brainf*ck code
# [-]>[-]<>+++++++[<++++++++++++++>-]<.++++++++++.+++..--.<

# [-]
while tape[pointer] != 0:
    tape[pointer] -= 1

# >
pointer += 1

# [-]
while tape[pointer] != 0:
    tape[pointer] -= 1

# <
pointer -= 1

# >
pointer += 1

# +++++++
tape[pointer] += 7

# [<++++++++++++++>-]
while tape[pointer] != 0:
    pointer -= 1
    tape[pointer] += 14
    pointer += 1
    tape[pointer] -= 1

# <
pointer -= 1

# .
output = chr(tape[pointer])

# ++++++++++
tape[pointer] += 10

# .
output += chr(tape[pointer])

# +++
tape[pointer] += 3

# .
output += chr(tape[pointer])

# .
output += chr(tape[pointer])

# --
tape[pointer] -= 2

# .
output += chr(tape[pointer])

# <
pointer -= 1

print(output)