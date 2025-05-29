def bf_interpreter(code):
    # Initialize tape with zeros and pointer
    tape = [0] * 30000
    ptr = 0
    output = ""
    i = 0
    
    # For debugging, let's also track the values
    debug_values = []
    
    while i < len(code):
        c = code[i]
        if c == '>':
            ptr += 1
        elif c == '<':
            ptr -= 1
        elif c == '+':
            tape[ptr] = (tape[ptr] + 1) % 256
        elif c == '-':
            tape[ptr] = (tape[ptr] - 1) % 256
        elif c == '.':
            output += chr(tape[ptr])
            debug_values.append(tape[ptr])  # Track ASCII values being output
        elif c == '[':
            if tape[ptr] == 0:
                brackets = 1
                while brackets:
                    i += 1
                    if code[i] == '[': brackets += 1
                    if code[i] == ']': brackets -= 1
        elif c == ']':
            if tape[ptr] != 0:
                brackets = 1
                while brackets:
                    i -= 1
                    if code[i] == '[': brackets -= 1
                    if code[i] == ']': brackets += 1
                i -= 1
        i += 1
    
    print(f"Output: {output}")
    print(f"ASCII values: {debug_values}")

code = "[-]>[-]<>++++++++[<+++++++++++++>-]<.---.----.>+++[<+++++++>-]<.>++++[<---->-]<-.<"
bf_interpreter(code)