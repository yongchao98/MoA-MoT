# Define the operations
def dollar(a, b):
    return (a + b) / b

def percent(a, b):
    return 4 * a + 4 * b

def hash(a, b):
    return (a + b) / b

# Calculate the expression step by step
step1 = dollar(3, 4)
step2 = percent(step1, 5)
result = hash(step2, 6)

print(result)