# Initialize memory tape and pointer
tape = [0] * 30000
pointer = 0

# Execute the Brainf*ck program
# >[-]>[-]<>+++++++[<++++++++++++++>-]<.-.>++++[<++++>-]<+.---.-.<

# >[-]
pointer += 1
tape[pointer] = 0

# >[-]
pointer += 1
tape[pointer] = 0

# <>
pointer -= 1

# +++++++
tape[pointer] += 7

# [<++++++++++++++>-]
while tape[pointer] != 0:
    pointer -= 1
    tape[pointer] += 14
    pointer += 1
    tape[pointer] -= 1

# <.-.
pointer -= 1
output = chr(tape[pointer])  # Output 'b'
tape[pointer] -= 1
output += chr(tape[pointer])  # Output 'a'

# >++++[<++++>-]
pointer += 1
tape[pointer] += 4
while tape[pointer] != 0:
    pointer -= 1
    tape[pointer] += 4
    pointer += 1
    tape[pointer] -= 1

# <+.
pointer -= 1
tape[pointer] += 1
output += chr(tape[pointer])  # Output 'r'

# ---.
tape[pointer] -= 3
output += chr(tape[pointer])  # Output 'o'

# -.
tape[pointer] -= 1
output += chr(tape[pointer])  # Output 'k'

# <.
pointer -= 1
output += chr(tape[pointer])  # Output 'e'

print(output)