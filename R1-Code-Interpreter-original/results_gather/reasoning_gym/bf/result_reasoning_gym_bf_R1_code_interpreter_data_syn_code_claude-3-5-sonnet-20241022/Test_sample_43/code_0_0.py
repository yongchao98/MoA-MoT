def bf_interpreter(code):
    # Initialize tape with zeros and pointer at position 0
    tape = [0] * 30000
    pointer = 0
    output = ""
    code_pointer = 0
    
    # Create bracket matching dictionary
    bracket_map = {}
    temp_stack = []
    for i, char in enumerate(code):
        if char == '[':
            temp_stack.append(i)
        elif char == ']':
            if temp_stack:
                start = temp_stack.pop()
                bracket_map[start] = i
                bracket_map[i] = start

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
                code_pointer = bracket_map[code_pointer]
        elif command == ']':
            if tape[pointer] != 0:
                code_pointer = bracket_map[code_pointer]
                
        code_pointer += 1
    
    return output

# The BF code to interpret
code = "[-]>[-]<>++++++++[<++++++++++++++>-]<.-.++++++.-.+++++.<"
result = bf_interpreter(code)
print(result)