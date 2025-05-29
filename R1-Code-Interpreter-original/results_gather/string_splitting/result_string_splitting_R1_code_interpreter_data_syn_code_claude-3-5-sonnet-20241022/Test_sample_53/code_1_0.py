def print_state(step, method, machines, parts):
    print(f"Step {step}, Method {method}:")
    print(f"Machines: A={machines['A']}, B={machines['B']}, C={machines['C']}")
    print(f"Parts: X={parts['X']}, Y={parts['Y']}, Z={parts['Z']}\n")

machines = {'A': 4, 'B': 2, 'C': 2}
parts = {'X': 0, 'Y': 0, 'Z': 0}
step = 0

# Detailed simulation with state printing
while True:
    changed = False
    for method in range(1, 6):
        if method == 1 and machines['A'] >= 1:
            machines['A'] -= 1
            parts['X'] += 2
            parts['Y'] += 1
            changed = True
            step += 1
            print_state(step, method, machines, parts)
        elif method == 2 and machines['B'] >= 2:
            machines['B'] -= 2
            parts['X'] += 1
            changed = True
            step += 1
            print_state(step, method, machines, parts)
        elif method == 3 and machines['C'] >= 2:
            machines['C'] -= 2
            parts['Y'] += 1
            changed = True
            step += 1
            print_state(step, method, machines, parts)
        elif method == 4 and machines['B'] >= 1 and machines['C'] >= 1:
            machines['B'] -= 1
            machines['C'] -= 1
            machines['A'] += 1
            changed = True
            step += 1
            print_state(step, method, machines, parts)
        elif method == 5 and parts['X'] >= 1 and parts['Y'] >= 1:
            parts['X'] -= 1
            parts['Y'] -= 1
            parts['Z'] += 1
            changed = True
            step += 1
            print_state(step, method, machines, parts)
    if not changed:
        break

result = [str(machines['A']), str(machines['B']), str(machines['C']), 
          str(parts['X']), str(parts['Y']), str(parts['Z'])]
print("Final result:", result)