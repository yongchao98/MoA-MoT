# Initial state
A, B, C = 1, 1, 2
X, Y, Z = 0, 0, 0

# Function to apply methods
def apply_methods(A, B, C, X, Y, Z):
    while True:
        # Track if any method was applied
        method_applied = False
        
        # Method 1
        if A >= 1:
            A -= 1
            X += 2
            Y += 1
            method_applied = True
        
        # Method 2
        elif B >= 2:
            B -= 2
            X += 1
            method_applied = True
        
        # Method 3
        elif C >= 2:
            C -= 2
            Y += 1
            method_applied = True
        
        # Method 4
        elif B >= 1 and C >= 1:
            B -= 1
            C -= 1
            A += 1
            method_applied = True
        
        # Method 5
        elif X >= 1 and Y >= 1:
            X -= 1
            Y -= 1
            Z += 1
            method_applied = True
        
        # If no method was applied, break the loop
        if not method_applied:
            break
    
    return A, B, C, X, Y, Z

# Apply methods
final_A, final_B, final_C, final_X, final_Y, final_Z = apply_methods(A, B, C, X, Y, Z)

# Output the final state
print([str(final_A), str(final_B), str(final_C), str(final_X), str(final_Y), str(final_Z)])