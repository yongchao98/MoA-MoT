# Initial counts
machines_A = 3
machines_B = 4
machines_C = 4
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to apply methods in sequence
def apply_methods():
    global machines_A, machines_B, machines_C, parts_X, parts_Y, parts_Z
    while True:
        # Track if any method was applied in this cycle
        method_applied = False
        
        # Method 1
        if machines_A >= 1:
            machines_A -= 1
            parts_X += 2
            parts_Y += 1
            method_applied = True
        
        # Method 2
        if machines_B >= 2:
            machines_B -= 2
            parts_X += 1
            method_applied = True
        
        # Method 3
        if machines_C >= 2:
            machines_C -= 2
            parts_Y += 1
            method_applied = True
        
        # Method 4
        if machines_B >= 1 and machines_C >= 1:
            machines_B -= 1
            machines_C -= 1
            machines_A += 1
            method_applied = True
        
        # Method 5
        if parts_X >= 1 and parts_Y >= 1:
            parts_X -= 1
            parts_Y -= 1
            parts_Z += 1
            method_applied = True
        
        # If no method was applied, break the loop
        if not method_applied:
            break

# Apply the methods
apply_methods()

# Output the final counts
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)