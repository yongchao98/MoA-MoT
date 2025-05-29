# Initialize the variables
A, B, C = 5, 5, 5
X = 0

# Define a flag to check if any method was executed in the cycle
method_executed = True

# Continue the cycle until no method can be executed
while method_executed:
    method_executed = False
    
    # Method 1: 1A + 1B + 1C -> 1X
    if A >= 1 and B >= 1 and C >= 1:
        A -= 1
        B -= 1
        C -= 1
        X += 1
        method_executed = True
    
    # Method 2: 1A + 1B -> 2C
    if A >= 1 and B >= 1:
        A -= 1
        B -= 1
        C += 2
        method_executed = True
    
    # Method 3: 2C -> 1X
    if C >= 2:
        C -= 2
        X += 1
        method_executed = True

# Output the final counts
print([A, B, C, X])