def brainfuck_interpreter(program):
    tape = [0] * 30000
    pointer = 0
    pc = 0
    output = []
    loop_stack = []

    while pc < len(program):
        command = program[pc]

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
        elif command == ',':
            pass  # No input handling
        elif command == '[':
            if tape[pointer] == 0:
                open_brackets = 1
                while open_brackets != 0:
                    pc += 1
                    if program[pc] == '[':
                        open_brackets += 1
                    elif program[pc] == ']':
                        open_brackets -= 1
            else:
                loop_stack.append(pc)
        elif command == ']':
            if tape[pointer] != 0:
                pc = loop_stack[-1]
            else:
                loop_stack.pop()

        pc += 1

    return ''.join(output)

program = "[-]>[-]<>++++++++[<+++++++++++++>-]<.+.+++++++++++.>++++[<---->-]<-.+++++.<"
print(brainfuck_interpreter(program))