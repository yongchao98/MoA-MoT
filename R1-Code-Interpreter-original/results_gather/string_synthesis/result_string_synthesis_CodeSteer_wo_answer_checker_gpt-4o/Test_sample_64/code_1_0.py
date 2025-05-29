def combine_blocks():
    # Initial counts of blocks
    count_A = 5
    count_B = 4
    count_C = 7
    count_A_curly = 0
    count_B_curly = 0
    count_C_curly = 0
    count_A_paren = 0
    count_B_paren = 0
    count_C_paren = 0

    # Function to apply each method
    def method_1():
        nonlocal count_A, count_B, count_C, count_A_curly
        if count_A >= 1 and count_B >= 1 and count_C >= 1:
            count_A -= 1
            count_B -= 1
            count_C -= 1
            count_A_curly += 1
            return True
        return False

    def method_2():
        nonlocal count_A, count_B, count_C_curly
        if count_A >= 1 and count_B >= 1:
            count_A -= 1
            count_B -= 1
            count_C_curly += 1
            return True
        return False

    def method_3():
        nonlocal count_B, count_C, count_B_curly
        if count_B >= 1 and count_C >= 1:
            count_B -= 1
            count_C -= 1
            count_B_curly += 1
            return True
        return False

    def method_4():
        nonlocal count_C, count_C_curly
        if count_C >= 2:
            count_C -= 2
            count_C_curly += 1
            return True
        return False

    def method_5():
        nonlocal count_A_curly, count_C_curly, count_A_paren, count_B_paren
        if count_A_curly >= 1 and count_C_curly >= 1:
            count_A_curly -= 1
            count_C_curly -= 1
            count_A_paren += 1
            count_B_paren += 1
            return True
        return False

    def method_6():
        nonlocal count_B_curly, count_C_paren
        if count_B_curly >= 2:
            count_B_curly -= 2
            count_C_paren += 1
            return True
        return False

    # Apply methods until no more changes can be made
    while True:
        changes_made = False
        for method in [method_1, method_2, method_3, method_4, method_5, method_6]:
            while method():
                changes_made = True
        if not changes_made:
            break

    # Output the final counts
    result = f"{count_A} [A], {count_B} [B], {count_C} [C], {count_A_curly} {{A}}, {count_B_curly} {{B}}, {count_C_curly} {{C}}, {count_A_paren} (A), {count_B_paren} (B), {count_C_paren} (C)"
    print(f"Your answer: {result}")

combine_blocks()