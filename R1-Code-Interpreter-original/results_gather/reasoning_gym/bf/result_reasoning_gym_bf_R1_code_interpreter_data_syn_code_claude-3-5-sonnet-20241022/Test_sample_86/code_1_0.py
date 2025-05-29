def bf_interpreter(code):
    # Initialize tape with zeros and pointer
    tape = [0] * 30000
    pointer = 0
    output = ""
    
    # Create bracket matching dictionary
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
    
    i = 0
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

# The BF code
code = "[-]>[-]<>+++++++++[<+++++++++++++>-]<+.-------------.--.++.+++."
bf_interpreter(code)