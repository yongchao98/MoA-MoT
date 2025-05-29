def bf_interpreter(code):
    tape = [0] * 30000  # Initialize tape with zeros
    pointer = 0         # Current position on tape
    output = ""         # Store output
    i = 0              # Code position
    
    # Create bracket matching dictionary for loops
    bracket_map = {}
    temp_stack = []
    for pos, char in enumerate(code):
        if char == '[':
            temp_stack.append(pos)
        elif char == ']':
            if temp_stack:
                start = temp_stack.pop()
                bracket_map[start] = pos
                bracket_map[pos] = start

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
                i = bracket_map[i]
        elif command == ']':
            if tape[pointer] != 0:
                i = bracket_map[i]
                
        i += 1
    
    print(output)

code = "[-]>[-]<>++++++++[<+++++++++++++>-]<-.------.>+++[<+++++++>-]<.>++++[<---->-]<-.+++++++.<"
bf_interpreter(code)