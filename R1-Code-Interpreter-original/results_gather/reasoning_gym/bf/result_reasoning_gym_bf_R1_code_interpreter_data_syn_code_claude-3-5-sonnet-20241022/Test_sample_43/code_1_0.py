def bf_interpreter(code):
    # Initialize tape with zeros and pointer at position 0
    tape = [0] * 30000
    pointer = 0
    output = ""
    
    # For debugging, let's also track the ASCII values
    ascii_values = []
    
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
            ascii_values.append(tape[pointer])
        elif command == '[':
            if tape[pointer] == 0:
                brackets = 1
                while brackets:
                    i += 1
                    if code[i] == '[': brackets += 1
                    if code[i] == ']': brackets -= 1
        elif command == ']':
            if tape[pointer] != 0:
                brackets = 1
                while brackets:
                    i -= 1
                    if code[i] == '[': brackets -= 1
                    if code[i] == ']': brackets += 1
                i -= 1
        i += 1
    
    print(f"Output: {output}")
    print(f"ASCII values: {ascii_values}")

code = "[-]>[-]<>++++++++[<++++++++++++++>-]<.-.++++++.-.+++++.<"
bf_interpreter(code)