# Initial counts
machines_A = 2
machines_B = 2
machines_C = 3
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to execute the methods in sequence
def execute_methods():
    global machines_A, machines_B, machines_C, parts_X, parts_Y, parts_Z
    while True:
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
        else:
            # If no method can be executed, break the loop
            break

# Execute the methods
execute_methods()

# Output the result
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)