# Initialize the quantities of machines and parts
machines_A = 4
machines_B = 2
machines_C = 2
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

# Loop through the methods until no method can be executed
while can_execute_any_method():
    # Method 1: One unit of machine A -> 2 units of part X and 1 unit of part Y
    if machines_A >= 1:
        machines_A -= 1
        parts_X += 2
        parts_Y += 1
    # Method 2: Two units of machine B -> 1 unit of part X
    elif machines_B >= 2:
        machines_B -= 2
        parts_X += 1
    # Method 3: Two machines C -> 1 part Y
    elif machines_C >= 2:
        machines_C -= 2
        parts_Y += 1
    # Method 4: One machine B and one machine C -> 1 machine A
    elif machines_B >= 1 and machines_C >= 1:
        machines_B -= 1
        machines_C -= 1
        machines_A += 1
    # Method 5: One part X and one part Y -> 1 part Z
    elif parts_X >= 1 and parts_Y >= 1:
        parts_X -= 1
        parts_Y -= 1
        parts_Z += 1

# Output the final quantities
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)