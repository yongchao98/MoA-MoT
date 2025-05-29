# Simulating the Brainf*ck program
def brainfuck_interpreter():
    tape = [0] * 30000
    pointer = 0
    code = "[-]>[-]<>+++++++++[<++++++++++++>-]<+.------------.+++++++++++++.-----.--------.<"
    output = []
    pc = 0

    while pc < len(code):
        command = code[pc]

        if command == '>':
            pointer += 1
        elif command == '<':
            pointer -= 1
        elif command == '+':
            tape[pointer] += 1
        elif command == '-':
            tape[pointer] -= 1
        elif command == '.':
            output.append(chr(tape[pointer]))
        elif command == '[':
            if tape[pointer] == 0:
                open_brackets = 1
                while open_brackets != 0:
                    pc += 1
                    if code[pc] == '[':
                        open_brackets += 1
                    elif code[pc] == ']':
                        open_brackets -= 1
        elif command == ']':
            if tape[pointer] != 0:
                close_brackets = 1
                while close_brackets != 0:
                    pc -= 1
                    if code[pc] == '[':
                        close_brackets -= 1
                    elif code[pc] == ']':
                        close_brackets += 1

        pc += 1

    return ''.join(output)

print(brainfuck_interpreter())