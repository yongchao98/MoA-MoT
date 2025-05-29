# Initialize the counts
machines_A = 3
machines_B = 3
machines_C = 3
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to check if any method can be executed
def can_execute_any_method():
    return (machines_A >= 1 or
            machines_B >= 2 or
            machines_C >= 2 or
            (machines_B >= 1 and machines_C >= 1) or
            (parts_X >= 1 and parts_Y >= 1))

# Cycle through the methods
while can_execute_any_method():
    # Method 1
    if machines_A >= 1:
        machines_A -= 1
        parts_X += 2
        parts_Y += 1
    # Method 2
    elif machines_B >= 2:
        machines_B -= 2
        parts_X += 1
    # Method 3
    elif machines_C >= 2:
        machines_C -= 2
        parts_Y += 1
    # Method 4
    elif machines_B >= 1 and machines_C >= 1:
        machines_B -= 1
        machines_C -= 1
        machines_A += 1
    # Method 5
    elif parts_X >= 1 and parts_Y >= 1:
        parts_X -= 1
        parts_Y -= 1
        parts_Z += 1

# Output the final counts
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)