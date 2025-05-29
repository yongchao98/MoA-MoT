def bf_interpret(code):
    memory = [0] * 30000
    pointer = 0
    output = ""
    i = 0
    
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
        elif command == '[':
            if memory[pointer] == 0:
                loop_count = 1
                while loop_count > 0:
                    i += 1
                    if code[i] == '[':
                        loop_count += 1
                    elif code[i] == ']':
                        loop_count -= 1
        elif command == ']':
            if memory[pointer] != 0:
                loop_count = 1
                while loop_count > 0:
                    i -= 1
                    if code[i] == '[':
                        loop_count -= 1
                    elif code[i] == ']':
                        loop_count += 1
                i -= 1
        i += 1
    
    print(output)

# The BF code to interpret
code = "[-]>[-]<>++++++++[<++++++++++++++>-]<++.-------------.+++++++++++++.+++.-------.<"
bf_interpret(code)