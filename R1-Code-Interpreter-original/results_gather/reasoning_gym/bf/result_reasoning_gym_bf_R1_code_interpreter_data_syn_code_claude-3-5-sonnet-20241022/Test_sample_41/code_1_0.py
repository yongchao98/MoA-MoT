def bf_interpreter(code):
    # Initialize memory with more precise tracking
    memory = [0] * 30000
    ptr = 0
    output = ""
    i = 0
    
    while i < len(code):
        c = code[i]
        if c == '>':
            ptr += 1
        elif c == '<':
            ptr -= 1
        elif c == '+':
            memory[ptr] = (memory[ptr] + 1) % 256
        elif c == '-':
            memory[ptr] = (memory[ptr] - 1) % 256
        elif c == '.':
            output += chr(memory[ptr])
        elif c == '[':
            if memory[ptr] == 0:
                brackets = 1
                while brackets:
                    i += 1
                    if code[i] == '[': brackets += 1
                    if code[i] == ']': brackets -= 1
        elif c == ']':
            if memory[ptr] != 0:
                brackets = 1
                while brackets:
                    i -= 1
                    if code[i] == '[': brackets -= 1
                    if code[i] == ']': brackets += 1
                i -= 1
        i += 1
    return output

code = "[-]>[-]<>+++++++++[<+++++++++++++>-]<--.++++.>+++[<------->-]<-.>++++[<++++>-]<+.-----.<"
print(bf_interpreter(code))