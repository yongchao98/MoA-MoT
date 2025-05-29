def bf_run(code):
    cells = [0] * 30000  # Memory cells
    cell_ptr = 0         # Cell pointer
    code_ptr = 0         # Code pointer
    output = []          # Store output characters
    
    # Find matching brackets
    brackets = {}
    stack = []
    for i, c in enumerate(code):
        if c == '[':
            stack.append(i)
        elif c == ']':
            if stack:
                start = stack.pop()
                brackets[start] = i
                brackets[i] = start
    
    # Execute code
    while code_ptr < len(code):
        c = code[code_ptr]
        
        if c == '>':
            cell_ptr += 1
        elif c == '<':
            cell_ptr -= 1
        elif c == '+':
            cells[cell_ptr] = (cells[cell_ptr] + 1) % 256
        elif c == '-':
            cells[cell_ptr] = (cells[cell_ptr] - 1) % 256
        elif c == '.':
            output.append(chr(cells[cell_ptr]))
        elif c == '[':
            if cells[cell_ptr] == 0:
                code_ptr = brackets[code_ptr]
        elif c == ']':
            if cells[cell_ptr] != 0:
                code_ptr = brackets[code_ptr]
        
        code_ptr += 1
    
    print(''.join(output))

# Run the code
code = "[-]>[-]<>++++++++[<++++++++++++++>-]<++.-------------.+++++++++++++.+++.-------.<"
bf_run(code)