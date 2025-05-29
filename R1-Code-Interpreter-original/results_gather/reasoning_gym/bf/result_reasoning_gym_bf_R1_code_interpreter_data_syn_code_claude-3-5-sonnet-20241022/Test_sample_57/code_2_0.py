def bf_interpreter(code):
    tape = [0] * 30000
    pointer = 0
    output = ""
    
    i = 0
    while i < len(code):
        cmd = code[i]
        
        if cmd == '>':
            pointer += 1
        elif cmd == '<':
            pointer -= 1
        elif cmd == '+':
            tape[pointer] = (tape[pointer] + 1) % 256
        elif cmd == '-':
            tape[pointer] = (tape[pointer] - 1) % 256
        elif cmd == '.':
            output += chr(tape[pointer])
            print(f"Output value: {tape[pointer]}")  # Debug print
        elif cmd == '[':
            if tape[pointer] == 0:
                brackets = 1
                while brackets:
                    i += 1
                    if code[i] == '[': brackets += 1
                    if code[i] == ']': brackets -= 1
        elif cmd == ']':
            if tape[pointer] != 0:
                brackets = 1
                while brackets:
                    i -= 1
                    if code[i] == '[': brackets -= 1
                    if code[i] == ']': brackets += 1
                continue
        i += 1
    
    print(f"Final tape values at output positions: {[tape[i] for i in range(5)]}")
    print(f"Output ASCII values: {[ord(c) for c in output]}")

code = "[-]>[-]<>++++++++++[<++++++++++>-]<++.-.+++++++++++++++.+.--.<"
bf_interpreter(code)