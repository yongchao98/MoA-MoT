def bf_interpreter(code):
    tape = [0] * 30000  # Initialize tape with zeros
    pointer = 0         # Current position on tape
    output = ""        # Store output
    code_pointer = 0   # Position in code
    
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
    
    print(output)

code = "[-]>[-]<>++++++++++[<++++++++++>-]<.+++++.+++++++++.++.+++++.<"
bf_interpreter(code)