def bf_analyze(code):
    # Initialize tape with zeros
    tape = [0] * 30000
    pointer = 0
    
    print("Step by step analysis:")
    # First part: [-]>[-]<> (clears first two cells and repositions)
    print("After clearing and repositioning, pointer at:", pointer)
    
    # Main calculation: +++++++++[<+++++++++++++>-]
    # 9 times outer loop, 13 times inner loop
    value = 9 * 13
    print("After multiplication loop, value:", value)
    
    # Character generation:
    # <--. (115 = 's')
    # -----. (110 = 'n')
    # +. (111 = 'o')
    # ++++++. (117 = 'u')
    # -. (116 = 't')
    
    result = ""
    base = value
    for adj in [-2, -5, 1, 6, -1]:
        result += chr(base + adj)
    print("Generated string:", result)

code = "[-]>[-]<>+++++++++[<+++++++++++++>-]<--.-----.+.++++++.-.<"
bf_analyze(code)