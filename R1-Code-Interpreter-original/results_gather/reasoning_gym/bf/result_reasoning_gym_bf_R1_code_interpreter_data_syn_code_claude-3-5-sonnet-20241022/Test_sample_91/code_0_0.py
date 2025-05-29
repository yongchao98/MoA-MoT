def bf_interpreter(code):
    # Initialize tape with zeros
    tape = [0] * 30000
    pointer = 0
    output = ""
    code_pointer = 0
    
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
            output += chr(tape[pointer])
        elif command == '[':
            if tape[pointer] == 0:
                bracket_count = 1
                while bracket_count > 0:
                    code_pointer += 1
                    if code[code_pointer] == '[':
                        bracket_count += 1
                    elif code[code_pointer] == ']':
                        bracket_count -= 1
        elif command == ']':
            if tape[pointer] != 0:
                bracket_count = 1
                while bracket_count > 0:
                    code_pointer -= 1
                    if code[code_pointer] == '[':
                        bracket_count -= 1
                    elif code[code_pointer] == ']':
                        bracket_count += 1
                continue
        
        code_pointer += 1
    
    return output

code = "[-]>[-]<>++++++++[<++++++++++++++>-]<.-.------.+++++.++++++.<"
print(repr(bf_interpreter(code)))