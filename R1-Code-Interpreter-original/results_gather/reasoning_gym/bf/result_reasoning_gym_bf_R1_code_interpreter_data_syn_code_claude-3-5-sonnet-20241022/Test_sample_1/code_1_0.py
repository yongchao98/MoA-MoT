def bf_interpreter(code):
    # Initialize memory with more detailed tracking
    memory = [0] * 30000
    pointer = 0
    output = ""
    
    def print_state():
        return f"Pointer: {pointer}, Current cell: {memory[pointer]}, ASCII: {chr(memory[pointer]) if 0 <= memory[pointer] <= 127 else '?'}"
    
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
    
    print(f"Final output: {repr(output)}")
    return output

code = "[-]>[-]<>++++++++++[<++++++++++>-]<+.-.+++++.---.>+++[<++++++>-]<+.<"
bf_interpreter(code)