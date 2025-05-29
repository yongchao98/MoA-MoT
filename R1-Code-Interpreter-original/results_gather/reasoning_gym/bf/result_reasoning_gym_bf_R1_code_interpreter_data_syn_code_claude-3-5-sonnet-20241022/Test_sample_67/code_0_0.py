def bf_interpreter(code):
    cells = [0] * 30000  # Initialize tape
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
            cells[pointer] = (cells[pointer] + 1) % 256
        elif command == '-':
            cells[pointer] = (cells[pointer] - 1) % 256
        elif command == '.':
            output += chr(cells[pointer])
        elif command == '[':
            if cells[pointer] == 0:
                bracket_count = 1
                while bracket_count > 0:
                    code_pointer += 1
                    if code[code_pointer] == '[':
                        bracket_count += 1
                    elif code[code_pointer] == ']':
                        bracket_count -= 1
        elif command == ']':
            if cells[pointer] != 0:
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

code = "[-]>[-]<>++++++++++[<+++++++++++>-]<+.--.--------.++.------.<"
bf_interpreter(code)