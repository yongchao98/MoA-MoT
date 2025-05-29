def bf_interpreter(code):
    # Initialize tape with zeros and debug info
    tape = [0] * 30000
    pointer = 0
    output = ""
    
    def print_debug(pos, val):
        return f"Position {pos}: {val} (ASCII: {chr(val) if 32 <= val <= 126 else '?'})"
    
    for i, cmd in enumerate(code):
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
        elif cmd == '[':
            if tape[pointer] == 0:
                loop_count = 1
                while loop_count > 0:
                    i += 1
                    if code[i] == '[': loop_count += 1
                    if code[i] == ']': loop_count -= 1
        elif cmd == ']':
            if tape[pointer] != 0:
                loop_count = 1
                while loop_count > 0:
                    i -= 1
                    if code[i] == '[': loop_count -= 1
                    if code[i] == ']': loop_count += 1
                i -= 1
    
    # Print the final values and output
    print(f"Final output: {output}")
    print(f"ASCII values: {[ord(c) for c in output]}")

code = "[-]>[-]<>++++++++++[<++++++++++>-]<++.-.+++++++++++++++.+.--.<"
bf_interpreter(code)