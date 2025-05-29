def bf_interpreter(code):
    tape = [0] * 30000
    pointer = 0
    output = ""
    
    # For debugging
    def print_state(msg):
        print(f"{msg}: cell[{pointer}] = {tape[pointer]}")
    
    i = 0
    while i < len(code):
        c = code[i]
        if c == '>':
            pointer += 1
        elif c == '<':
            pointer -= 1
        elif c == '+':
            tape[pointer] = (tape[pointer] + 1) % 256
        elif c == '-':
            tape[pointer] = (tape[pointer] - 1) % 256
        elif c == '.':
            output += chr(tape[pointer])
            print(f"Output char: {chr(tape[pointer])} (ASCII: {tape[pointer]})")
        elif c == '[':
            if tape[pointer] == 0:
                brackets = 1
                while brackets:
                    i += 1
                    if code[i] == '[': brackets += 1
                    if code[i] == ']': brackets -= 1
        elif c == ']':
            if tape[pointer] != 0:
                brackets = 1
                while brackets:
                    i -= 1
                    if code[i] == '[': brackets -= 1
                    if code[i] == ']': brackets += 1
                i -= 1
        i += 1
    return output

code = "[-]>[-]<>+++++++++[<+++++++++++++>-]<--.-----.+.++++++.-.<"
result = bf_interpreter(code)
print(f"\nFinal output: {result}")