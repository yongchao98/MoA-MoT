# Initialize the quantities of machines and parts
A, B, C = 1, 1, 2
X, Y, Z = 0, 0, 0

# Function to simulate the process
def simulate_dismantling():
    global A, B, C, X, Y, Z
    while True:
        # Track if any method was executed in this cycle
        executed = False
        
        # Method 1
        if A >= 1:
            A -= 1
            X += 2
            Y += 1
            executed = True
        
        # Method 2
        elif B >= 2:
            B -= 2
            X += 1
            executed = True
        
        # Method 3
        elif C >= 2:
            C -= 2
            Y += 1
            executed = True
        
        # Method 4
        elif B >= 1 and C >= 1:
            B -= 1
            C -= 1
            A += 1
            executed = True
        
        # Method 5
        elif X >= 1 and Y >= 1:
            X -= 1
            Y -= 1
            Z += 1
            executed = True
        
        # If no method was executed, break the loop
        if not executed:
            break

# Run the simulation
simulate_dismantling()

# Output the final quantities
result = [str(A), str(B), str(C), str(X), str(Y), str(Z)]
print(result)