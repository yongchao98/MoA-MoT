def bf_interpreter(code):
    # Initialize memory array and pointers
    memory = [0] * 30000
    pointer = 0
    code_pointer = 0
    output = ""
    
    # For debugging, let's also track the ASCII values
    ascii_values = []
    
    while code_pointer < len(code):
        command = code[code_pointer]
        
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
            ascii_values.append(memory[pointer])  # Track ASCII values
        elif command == '[':
            if memory[pointer] == 0:
                bracket_count = 1
                while bracket_count > 0:
                    code_pointer += 1
                    if code[code_pointer] == '[': bracket_count += 1
                    if code[code_pointer] == ']': bracket_count -= 1
        elif command == ']':
            if memory[pointer] != 0:
                bracket_count = 1
                while bracket_count > 0:
                    code_pointer -= 1
                    if code[code_pointer] == '[': bracket_count -= 1
                    if code[code_pointer] == ']': bracket_count += 1
                code_pointer -= 1
        
        code_pointer += 1
    
    print(f"ASCII values of output: {ascii_values}")
    return output

code = "[-]>[-]<>++++++++[<++++++++++++>-]<+.>++++[<++++++>-]<-.---------------.++++++.-.<<"
print(bf_interpreter(code))