# Initialize the memory tape and data pointer
tape = [0] * 30000
ptr = 0

# Execute the Brainf*ck code
code = "[-]>[-]<>++++++++++[<+++++++++++>-]<+.+++++++++.---------------.-----.+.<"
i = 0
output = []

while i < len(code):
    cmd = code[i]
    if cmd == '>':
        ptr += 1
    elif cmd == '<':
        ptr -= 1
    elif cmd == '+':
        tape[ptr] += 1
    elif cmd == '-':
        tape[ptr] -= 1
    elif cmd == '.':
        output.append(chr(tape[ptr]))
    elif cmd == '[':
        if tape[ptr] == 0:
            open_brackets = 1
            while open_brackets != 0:
                i += 1
                if code[i] == '[':
                    open_brackets += 1
                elif code[i] == ']':
                    open_brackets -= 1
    elif cmd == ']':
        if tape[ptr] != 0:
            close_brackets = 1
            while close_brackets != 0:
                i -= 1
                if code[i] == '[':
                    close_brackets -= 1
                elif code[i] == ']':
                    close_brackets += 1
    i += 1

print(''.join(output))