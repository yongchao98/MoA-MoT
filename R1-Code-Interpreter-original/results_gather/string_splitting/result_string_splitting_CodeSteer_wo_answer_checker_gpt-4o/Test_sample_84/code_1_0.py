# Initialize the counts of machines and parts
machines_A = 5
machines_B = 3
machines_C = 6
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to execute the methods in sequence
def execute_methods():
    global machines_A, machines_B, machines_C, parts_X, parts_Y, parts_Z
    while True:
        # Method 1: Dismantle 1 machine A -> 2 parts X, 1 part Y
        if machines_A >= 1:
            machines_A -= 1
            parts_X += 2
            parts_Y += 1
        # Method 2: Dismantle 2 machines B -> 1 part X
        elif machines_B >= 2:
            machines_B -= 2
            parts_X += 1
        # Method 3: Dismantle 2 machines C -> 1 part Y
        elif machines_C >= 2:
            machines_C -= 2
            parts_Y += 1
        # Method 4: Combine 1 machine B and 1 machine C -> 1 machine A
        elif machines_B >= 1 and machines_C >= 1:
            machines_B -= 1
            machines_C -= 1
            machines_A += 1
        # Method 5: Combine 1 part X and 1 part Y -> 1 part Z
        elif parts_X >= 1 and parts_Y >= 1:
            parts_X -= 1
            parts_Y -= 1
            parts_Z += 1
        else:
            # No more methods can be executed
            break

# Execute the methods
execute_methods()

# Output the final counts
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)