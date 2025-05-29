def can_execute_method(method, machines, parts):
    if method == 1:
        return machines[0] >= 1  # Need 1 machine A
    elif method == 2:
        return machines[1] >= 2  # Need 2 machine B
    elif method == 3:
        return machines[2] >= 2  # Need 2 machine C
    elif method == 4:
        return machines[1] >= 1 and machines[2] >= 1  # Need 1 B and 1 C
    elif method == 5:
        return parts[0] >= 1 and parts[1] >= 1  # Need 1 X and 1 Y
    return False

def execute_method(method, machines, parts):
    if method == 1:
        machines[0] -= 1
        parts[0] += 2  # Get 2 X
        parts[1] += 1  # Get 1 Y
    elif method == 2:
        machines[1] -= 2
        parts[0] += 1  # Get 1 X
    elif method == 3:
        machines[2] -= 2
        parts[1] += 1  # Get 1 Y
    elif method == 4:
        machines[1] -= 1
        machines[2] -= 1
        machines[0] += 1  # Get 1 A
    elif method == 5:
        parts[0] -= 1
        parts[1] -= 1
        parts[2] += 1  # Get 1 Z

def solve_dismantling():
    # Initialize [A, B, C] machines and [X, Y, Z] parts
    machines = [2, 2, 4]
    parts = [0, 0, 0]
    current_method = 1
    
    while True:
        # Try to execute current method
        if can_execute_method(current_method, machines, parts):
            execute_method(current_method, machines, parts)
        else:
            # Try next method
            current_method = current_method % 5 + 1
            # Check if we can execute any method
            can_continue = False
            for m in range(1, 6):
                if can_execute_method(m, machines, parts):
                    can_continue = True
                    break
            if not can_continue:
                break
    
    # Convert results to strings
    result = [str(x) for x in machines + parts]
    print(result)

solve_dismantling()