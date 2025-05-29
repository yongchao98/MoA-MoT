def bf(code):
    mem = [0] * 30000
    p = 0
    i = 0
    out = []
    while i < len(code):
        c = code[i]
        if c == '>': p += 1
        elif c == '<': p -= 1
        elif c == '+': mem[p] = (mem[p] + 1) % 256
        elif c == '-': mem[p] = (mem[p] - 1) % 256
        elif c == '.': out.append(chr(mem[p]))
        elif c == '[':
            if mem[p] == 0:
                nest = 1
                while nest:
                    i += 1
                    if code[i] == '[': nest += 1
                    if code[i] == ']': nest -= 1
        elif c == ']':
            if mem[p]:
                nest = 1
                while nest:
                    i -= 1
                    if code[i] == '[': nest -= 1
                    if code[i] == ']': nest += 1
                continue
        i += 1
    print(''.join(out), end='')

bf("[-]>[-]<>++++++++[<++++++++++++++>-]<++.-------------.+++++++++++++.+++.-------.<")