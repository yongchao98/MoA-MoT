def bf_interpreter(code):
    # Initialize memory with 30000 cells
    memory = [0] * 30000
    pointer = 0
    output = ""
    i = 0
    
    # For debugging, let's also track cell values
    debug_cells = []
    
    while i < len(code):
        command = code[i]
        
        if command == '>':
            pointer += 1
        elif command == '<':
            pointer -= 1
        elif command == '+':
            memory[pointer] = (memory[pointer] + 1) % 256
        elif command == '-':
            memory[pointer] = (memory[pointer] - 1) % 256
        elif command == '.':
            output += chr(memory[pointer])
            debug_cells.append(memory[pointer])  # Track output values
        elif command == '[':
            if memory[pointer] == 0:
                bracket_count = 1
                while bracket_count > 0:
                    i += 1
                    if code[i] == '[': bracket_count += 1
                    if code[i] == ']': bracket_count -= 1
        elif command == ']':
            if memory[pointer] != 0:
                bracket_count = 1
                while bracket_count > 0:
                    i -= 1
                    if code[i] == '[': bracket_count -= 1
                    if code[i] == ']': bracket_count += 1
                i -= 1
        
        i += 1
    
    print(f"Output: {output}")
    print(f"ASCII values of output: {debug_cells}")

program = "[-]>[-]<>+++++++++[<+++++++++++++>-]<+.>++++[<---->-]<-.+++++++++.+++++++.>++++[<---->-]<.<"
bf_interpreter(program)