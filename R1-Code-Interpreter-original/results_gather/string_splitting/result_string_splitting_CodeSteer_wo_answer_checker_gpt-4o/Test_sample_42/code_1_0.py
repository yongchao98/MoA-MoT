# Initial quantities
machines_A = 2
machines_B = 2
machines_C = 4
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to execute the methods in order
def execute_methods():
    global machines_A, machines_B, machines_C, parts_X, parts_Y, parts_Z
    while True:
        executed = False
        
        # Method 1: Dismantle machine A
        if machines_A >= 1:
            machines_A -= 1
            parts_X += 2
            parts_Y += 1
            executed = True
        
        # Method 2: Dismantle machine B
        elif machines_B >= 2:
            machines_B -= 2
            parts_X += 1
            executed = True
        
        # Method 3: Dismantle machine C
        elif machines_C >= 2:
            machines_C -= 2
            parts_Y += 1
            executed = True
        
        # Method 4: Combine machine B and C into machine A
        elif machines_B >= 1 and machines_C >= 1:
            machines_B -= 1
            machines_C -= 1
            machines_A += 1
            executed = True
        
        # Method 5: Combine parts X and Y into part Z
        elif parts_X >= 1 and parts_Y >= 1:
            parts_X -= 1
            parts_Y -= 1
            parts_Z += 1
            executed = True
        
        # If no method was executed, break the loop
        if not executed:
            break

# Execute the methods
execute_methods()

# Output the final quantities
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)