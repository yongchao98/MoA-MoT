def bf_interpreter(code):
    tape = [0] * 30000  # Initialize tape with zeros
    pointer = 0         # Data pointer
    output = ""        # Store output
    i = 0              # Code pointer
    
    # Create bracket matching dictionary
    brackets = {}
    stack = []
    for pos, char in enumerate(code):
        if char == '[':
            stack.append(pos)
        elif char == ']':
            if stack:
                start = stack.pop()
                brackets[start] = pos
                brackets[pos] = start

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
        elif command == '[':
            if tape[pointer] == 0:
                i = brackets[i]
        elif command == ']':
            if tape[pointer] != 0:
                i = brackets[i]
        
        i += 1
    
    print(output)

# The BF code to interpret
code = "[-]>[-]<>++++++++[<+++++++++++++>-]<-.+++++++++++.---.++++++++.-----------.<"
bf_interpreter(code)