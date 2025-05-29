def bf_interpreter(code):
    # Initialize memory and pointer
    memory = [0] * 30000
    pointer = 0
    output = ""
    
    # For debugging, let's also track the ASCII values
    ascii_values = []
    
    # Process each character
    i = 0
    while i < len(code):
        c = code[i]
        if c == '>':
            pointer += 1
        elif c == '<':
            pointer -= 1
        elif c == '+':
            memory[pointer] = (memory[pointer] + 1) % 256
        elif c == '-':
            memory[pointer] = (memory[pointer] - 1) % 256
        elif c == '.':
            output += chr(memory[pointer])
            ascii_values.append(memory[pointer])
        elif c == '[':
            if memory[pointer] == 0:
                brackets = 1
                while brackets:
                    i += 1
                    if code[i] == '[': brackets += 1
                    if code[i] == ']': brackets -= 1
        elif c == ']':
            if memory[pointer] != 0:
                brackets = 1
                while brackets:
                    i -= 1
                    if code[i] == '[': brackets -= 1
                    if code[i] == ']': brackets += 1
                i -= 1
        i += 1
    
    print(f"ASCII values: {ascii_values}")
    return output

code = "[-]>[-]<>+++++++++[<+++++++++++++>-]<+.-------------.+++++++++.>++++[<---->-]<-.+++++++++++.<"
print(bf_interpreter(code))