# Initial counts
machines_A = 5
machines_B = 5
machines_C = 4
parts_X = 0
parts_Y = 0
parts_Z = 0

# Function to simulate the process
def simulate_dismantling():
    global machines_A, machines_B, machines_C, parts_X, parts_Y, parts_Z
    while True:
        # Method 1: Dismantle a single machine A
        if machines_A >= 1:
            machines_A -= 1
            parts_X += 2
            parts_Y += 1
        # Method 2: Dismantle two machines B
        elif machines_B >= 2:
            machines_B -= 2
            parts_X += 1
        # Method 3: Dismantle two machines C
        elif machines_C >= 2:
            machines_C -= 2
            parts_Y += 1
        # Method 4: Combine one machine B and one machine C
        elif machines_B >= 1 and machines_C >= 1:
            machines_B -= 1
            machines_C -= 1
            machines_A += 1
        # Method 5: Combine one part X and one part Y
        elif parts_X >= 1 and parts_Y >= 1:
            parts_X -= 1
            parts_Y -= 1
            parts_Z += 1
        else:
            # No more methods can be executed
            break

# Run the simulation
simulate_dismantling()

# Output the final counts
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)