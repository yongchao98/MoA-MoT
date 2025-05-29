# Initial quantities
machines_A = 3
machines_B = 2
machines_C = 1
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to execute the methods in order
def execute_methods():
    global machines_A, machines_B, machines_C, parts_X, parts_Y, parts_Z
    while True:
        # Method 1: One unit of machine A can be dismantled into two units of part X and one unit of part Y
        if machines_A >= 1:
            machines_A -= 1
            parts_X += 2
            parts_Y += 1
            continue
        
        # Method 2: Every round, dismantling two machines B results in one part X
        if machines_B >= 2:
            machines_B -= 2
            parts_X += 1
            continue
        
        # Method 3: Each round, two machines C break down into one part Y
        if machines_C >= 2:
            machines_C -= 2
            parts_Y += 1
            continue
        
        # Method 4: Every round, merging one machine B with one machine C creates one machine A
        if machines_B >= 1 and machines_C >= 1:
            machines_B -= 1
            machines_C -= 1
            machines_A += 1
            continue
        
        # Method 5: Each round, one part X and one part Y merge to form one part Z
        if parts_X >= 1 and parts_Y >= 1:
            parts_X -= 1
            parts_Y -= 1
            parts_Z += 1
            continue
        
        # If no method can be executed, break the loop
        break

# Execute the methods
execute_methods()

# Output the final quantities
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)