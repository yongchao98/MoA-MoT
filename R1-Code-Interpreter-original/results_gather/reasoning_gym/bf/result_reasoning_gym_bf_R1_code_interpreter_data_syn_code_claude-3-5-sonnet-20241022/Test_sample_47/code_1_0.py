def bf_interpreter(code):
    # Initialize memory and pointers
    memory = [0] * 30000
    ptr = 0
    output = ""
    i = 0
    
    # For debugging
    def print_relevant_memory():
        return [memory[j] for j in range(max(0, ptr-2), ptr+3)]
    
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
            print(f"Output char: {chr(memory[ptr])} (ASCII: {memory[ptr]})")
        elif c == '[':
            if memory[ptr] == 0:
                count = 1
                while count > 0:
                    i += 1
                    if code[i] == '[': count += 1
                    if code[i] == ']': count -= 1
        elif c == ']':
            if memory[ptr] != 0:
                count = 1
                while count > 0:
                    i -= 1
                    if code[i] == '[': count -= 1
                    if code[i] == ']': count += 1
                i -= 1
        i += 1
    print(output)

code = "[-]>[-]<>++++++++[<+++++++++++++>-]<.+++++++++++++.--------.------------.+++++++++++++.<"
bf_interpreter(code)