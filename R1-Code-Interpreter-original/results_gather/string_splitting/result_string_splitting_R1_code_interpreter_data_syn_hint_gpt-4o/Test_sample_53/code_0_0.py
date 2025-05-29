# Initial counts
machines_A = 4
machines_B = 2
machines_C = 2
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to simulate the process
def simulate_dismantling():
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
            break

simulate_dismantling()

# Output the final counts
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)