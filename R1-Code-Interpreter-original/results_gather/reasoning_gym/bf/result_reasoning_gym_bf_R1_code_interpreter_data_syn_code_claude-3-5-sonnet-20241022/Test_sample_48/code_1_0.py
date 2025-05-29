def analyze_bf_values(code):
    # Let's track the cell values at each step
    tape = [0] * 30000
    pointer = 0
    steps = []
    
    # First part: [-]>[-]<> (clears first two cells and repositions)
    # Then: +++++++++[<+++++++++++++>-]
    # This creates: 9 loops of adding 13 to the first cell = 117
    # After that: <--. (115 = 's')
    # Then: -----------. (104 = 'h')
    # Then: ---. (101 = 'e')
    # Then: ++++. (105 = 'i')
    # Finally: ++. (107 = 'k')
    
    base = 9 * 13  # Initial calculation
    print(f"Base calculation: {base}")
    print(f"After first decrease: {base-2} ('{chr(base-2)}')")
    print(f"After second set: {base-2-11} ('{chr(base-2-11)}')")
    print(f"After third set: {base-2-11-3} ('{chr(base-2-11-3)}')")
    print(f"After fourth set: {base-2-11-3+4} ('{chr(base-2-11-3+4)}')")
    print(f"After final set: {base-2-11-3+4+2} ('{chr(base-2-11-3+4+2)}')")

analyze_bf_values("[-]>[-]<>+++++++++[<+++++++++++++>-]<--.-----------.---.++++.++.<")