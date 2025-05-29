def bf_interpreter(code):
    # Initialize tape with zeros and pointer at position 0
    tape = [0] * 30000
    pointer = 0
    output = ""
    i = 0
    
    # Create bracket matching dictionary
    bracket_map = {}
    stack = []
    for pos, char in enumerate(code):
        if char == '[':
            stack.append(pos)
        elif char == ']':
            if stack:
                start = stack.pop()
                bracket_map[start] = pos
                bracket_map[pos] = start
    
    # For debugging
    while i < len(code):
        command = code[i]
        
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
            # print(f"Output char: {chr(tape[pointer])} (ASCII: {tape[pointer]})")  # Debug
        elif command == '[':
            if tape[pointer] == 0:
                i = bracket_map[i]
        elif command == ']':
            if tape[pointer] != 0:
                i = bracket_map[i]
        
        i += 1
    
    print(output)

code = "[-]>[-]<>+++++++++[<+++++++++++>-]<.>+++[<++++++>-]<.---.++++.>++++[<---->-]<-.<"
bf_interpreter(code)