def brainfuck_interpreter(code):
    tape = [0] * 30000
    pointer = 0
    code_pointer = 0
    output = []
    loop_stack = []

    while code_pointer < len(code):
        command = code[code_pointer]

        if command == '>':
            pointer += 1
        elif command == '<':
            pointer -= 1
        elif command == '+':
            tape[pointer] = (tape[pointer] + 1) % 256
        elif command == '-':
            tape[pointer] = (tape[pointer] - 1) % 256
        elif command == '.':
            output.append(chr(tape[pointer]))
        elif command == '[':
            if tape[pointer] == 0:
                open_brackets = 1
                while open_brackets != 0:
                    code_pointer += 1
                    if code[code_pointer] == '[':
                        open_brackets += 1
                    elif code[code_pointer] == ']':
                        open_brackets -= 1
            else:
                loop_stack.append(code_pointer)
        elif command == ']':
            if tape[pointer] != 0:
                code_pointer = loop_stack[-1]
            else:
                loop_stack.pop()

        code_pointer += 1

    return ''.join(output)

bf_code = "[-]>[-]<>+++++++++[<+++++++++++++>-]<+.-------------.+++++++++.>++++[<---->-]<-.+++++++++++.<"
print(brainfuck_interpreter(bf_code))