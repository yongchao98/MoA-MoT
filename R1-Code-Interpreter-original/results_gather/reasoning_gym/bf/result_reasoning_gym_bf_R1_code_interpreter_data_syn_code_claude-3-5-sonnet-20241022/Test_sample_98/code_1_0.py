def bf_interpreter(code):
    # Initialize tape and debug info
    tape = [0] * 30000
    pointer = 0
    output = ""
    code_pointer = 0
    
    # For debugging
    states = []
    
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
            # Store state at each output
            states.append(f"Output '{chr(tape[pointer])}' (ASCII {tape[pointer]})")
        elif command == '[':
            if tape[pointer] == 0:
                code_pointer = brackets[code_pointer]
        elif command == ']':
            if tape[pointer] != 0:
                code_pointer = brackets[code_pointer]
        
        code_pointer += 1
    
    # Print debug info
    for state in states:
        print(state)
    print("\nFinal output:", output)

code = "[-]>[-]<>++++++++[<+++++++++++++>-]<.+++++++++++++.---..+++++++.<"
bf_interpreter(code)