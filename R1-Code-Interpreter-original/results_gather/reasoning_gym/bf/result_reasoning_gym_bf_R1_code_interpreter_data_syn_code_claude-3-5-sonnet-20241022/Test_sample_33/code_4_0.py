def run_bf(code):
    cells = [0] * 30000
    cell_ptr = 0
    i = 0
    output = ""
    
    while i < len(code):
        c = code[i]
        if c == '>':
            cell_ptr += 1
        elif c == '<':
            cell_ptr -= 1
        elif c == '+':
            cells[cell_ptr] = (cells[cell_ptr] + 1) % 256
        elif c == '-':
            cells[cell_ptr] = (cells[cell_ptr] - 1) % 256
        elif c == '.':
            output += chr(cells[cell_ptr])
        elif c == '[' and cells[cell_ptr] == 0:
            count = 1
            while count > 0:
                i += 1
                if code[i] == '[': count += 1
                if code[i] == ']': count -= 1
        elif c == ']' and cells[cell_ptr] != 0:
            count = 1
            while count > 0:
                i -= 1
                if code[i] == '[': count -= 1
                if code[i] == ']': count += 1
        i += 1
    print(output, end='')

code = "[-]>[-]<>++++++++[<++++++++++++++>-]<++.-------------.+++++++++++++.+++.-------.<"
run_bf(code)