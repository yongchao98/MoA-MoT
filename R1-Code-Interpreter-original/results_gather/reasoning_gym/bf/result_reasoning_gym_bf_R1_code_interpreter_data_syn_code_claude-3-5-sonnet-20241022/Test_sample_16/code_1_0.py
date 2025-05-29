def bf_interpreter(code):
    memory = [0] * 30000
    pointer = 0
    output = ""
    code_pos = 0
    
    while code_pos < len(code):
        command = code[code_pos]
        
        if command == '>':
            pointer += 1
        elif command == '<':
            pointer -= 1
        elif command == '+':
            memory[pointer] = (memory[pointer] + 1) % 256
        elif command == '-':
            memory[pointer] = (memory[pointer] - 1) % 256
        elif command == '.':
            output += chr(memory[pointer])
        elif command == '[':
            if memory[pointer] == 0:
                bracket_count = 1
                while bracket_count > 0:
                    code_pos += 1
                    if code[code_pos] == '[':
                        bracket_count += 1
                    elif code[code_pos] == ']':
                        bracket_count -= 1
            else:
                pass
        elif command == ']':
            if memory[pointer] != 0:
                bracket_count = 1
                while bracket_count > 0:
                    code_pos -= 1
                    if code[code_pos] == '[':
                        bracket_count -= 1
                    elif code[code_pos] == ']':
                        bracket_count += 1
                code_pos -= 1
        
        code_pos += 1
    
    print(output)

code = "[-]>[-]<>+++++++++[<++++++++++++>-]<+.++.++++..++++++."
bf_interpreter(code)