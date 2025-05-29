def bf_interpreter_debug(code):
    tape = [0] * 30000
    pointer = 0
    output = ""
    code_pointer = 0
    
    # Find matching brackets
    brackets = {}
    stack = []
    for i, char in enumerate(code):
        if char == '[':
            stack.append(i)
        elif char == ']':
            if stack:
                start = stack.pop()
                brackets[start] = i
                brackets[i] = start
    
    # For debugging
    outputs_with_ascii = []
    
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
            outputs_with_ascii.append((chr(tape[pointer]), tape[pointer]))
            output += chr(tape[pointer])
        elif command == '[':
            if tape[pointer] == 0:
                code_pointer = brackets[code_pointer]
        elif command == ']':
            if tape[pointer] != 0:
                code_pointer = brackets[code_pointer]
        
        code_pointer += 1
    
    print("Characters and their ASCII values:")
    for char, ascii_val in outputs_with_ascii:
        print(f"'{char}' (ASCII: {ascii_val})")
    print("\nFinal output:", output)

code = "[-]>[-]<>+++++++++[<+++++++++++++>-]<--.-----------.+++++++.+++.----.<"
bf_interpreter_debug(code)