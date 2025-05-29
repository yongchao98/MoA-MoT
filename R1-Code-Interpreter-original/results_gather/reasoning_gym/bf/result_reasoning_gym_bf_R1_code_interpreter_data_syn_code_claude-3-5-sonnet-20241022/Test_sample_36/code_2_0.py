def bf_interpreter(code):
    memory = [0] * 30000
    pointer = 0
    output = ""
    
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
            print(f"Output character: {chr(memory[pointer])} (ASCII: {memory[pointer]})")
        elif c == '[':
            if memory[pointer] == 0:
                loop_count = 1
                while loop_count > 0:
                    i += 1
                    if code[i] == '[': loop_count += 1
                    if code[i] == ']': loop_count -= 1
        elif c == ']':
            if memory[pointer] != 0:
                loop_count = 1
                while loop_count > 0:
                    i -= 1
                    if code[i] == '[': loop_count -= 1
                    if code[i] == ']': loop_count += 1
                continue
        
        i += 1
        if c in '+-><':
            print(f"Memory[{pointer}] = {memory[pointer]}")
    
    return output

code = "[-]>[-]<>+++++++++[<+++++++++++++>-]<--.+.-----.-.---------.<"
result = bf_interpreter(code)
print("\nFinal output:", repr(result))