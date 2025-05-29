def bf_interpreter(code):
    # Initialize tape and pointer
    tape = [0] * 30000
    ptr = 0
    output = ""
    
    # For debugging, track cell values
    def print_cells():
        return f"Cells[0:5]: {tape[0:5]}, Pointer at: {ptr}, Value: {tape[ptr]}"
    
    i = 0
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
            print(f"Output char: {chr(tape[ptr])} (ASCII: {tape[ptr]})")
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
        
    print(f"\nFinal output: {output}")

code = "[-]>[-]<>++++++++++[<+++++++++++>-]<+.+.-------.++++++++++++.--------.<"
bf_interpreter(code)