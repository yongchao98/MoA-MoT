# Initial counts
machines_A = 2
machines_B = 4
machines_C = 3
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to execute the methods in sequence
def execute_methods():
    global machines_A, machines_B, machines_C, parts_X, parts_Y, parts_Z
    executed = False

    # Method 1
    if machines_A >= 1:
        machines_A -= 1
        parts_X += 2
        parts_Y += 1
        executed = True

    # Method 2
    if machines_B >= 2:
        machines_B -= 2
        parts_X += 1
        executed = True

    # Method 3
    if machines_C >= 2:
        machines_C -= 2
        parts_Y += 1
        executed = True

    # Method 4
    if machines_B >= 1 and machines_C >= 1:
        machines_B -= 1
        machines_C -= 1
        machines_A += 1
        executed = True

    # Method 5
    if parts_X >= 1 and parts_Y >= 1:
        parts_X -= 1
        parts_Y -= 1
        parts_Z += 1
        executed = True

    return executed

# Loop until no methods can be executed
while execute_methods():
    pass

# Output the final counts
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)