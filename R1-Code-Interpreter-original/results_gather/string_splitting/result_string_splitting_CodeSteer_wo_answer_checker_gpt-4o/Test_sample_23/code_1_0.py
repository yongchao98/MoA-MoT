# Initialize the counts
machines_A = 3
machines_B = 2
machines_C = 2
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
            continue
        
        # Method 2
        if machines_B >= 2:
            machines_B -= 2
            parts_X += 1
            continue
        
        # Method 3
        if machines_C >= 2:
            machines_C -= 2
            parts_Y += 1
            continue
        
        # Method 4
        if machines_B >= 1 and machines_C >= 1:
            machines_B -= 1
            machines_C -= 1
            machines_A += 1
            continue
        
        # Method 5
        if parts_X >= 1 and parts_Y >= 1:
            parts_X -= 1
            parts_Y -= 1
            parts_Z += 1
            continue
        
        # If no method can be executed, break the loop
        break

# Execute the methods
execute_methods()

# Output the final counts
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)