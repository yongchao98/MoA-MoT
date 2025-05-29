def bf_interpreter(code):
    # Initialize tape with zeros
    tape = [0] * 30000
    pointer = 0
    output = ""
    code_pointer = 0
    
    # For debugging
    debug_values = []
    
    # Create bracket matching dictionary
    bracket_map = {}
    stack = []
    for i, char in enumerate(code):
        if char == '[':
            stack.append(i)
        elif char == ']':
            start = stack.pop()
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
            debug_values.append(tape[pointer])  # Store ASCII values for debugging
        elif command == '[':
            if tape[pointer] == 0:
                code_pointer = bracket_map[code_pointer]
        elif command == ']':
            if tape[pointer] != 0:
                code_pointer = bracket_map[code_pointer]
        
        code_pointer += 1
    
    print(f"Debug ASCII values: {debug_values}")  # Print ASCII values for verification
    return output

code = "[-]>[-]<>+++++++[<+++++++++++++++>-]<.-----.+++++.++++++.--.<"
print(repr(bf_interpreter(code)))