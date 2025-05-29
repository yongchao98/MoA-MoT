# Initial quantities
machines_A = 4
machines_B = 3
machines_C = 5
parts_X = 0
parts_Y = 0
parts_Z = 0

# Define the methods
def method_1():
    global machines_A, parts_X, parts_Y
    if machines_A >= 1:
        machines_A -= 1
        parts_X += 2
        parts_Y += 1
        return True
    return False

def method_2():
    global machines_B, parts_X
    if machines_B >= 2:
        machines_B -= 2
        parts_X += 1
        return True
    return False

def method_3():
    global machines_C, parts_Y
    if machines_C >= 2:
        machines_C -= 2
        parts_Y += 1
        return True
    return False

def method_4():
    global machines_B, machines_C, machines_A
    if machines_B >= 1 and machines_C >= 1:
        machines_B -= 1
        machines_C -= 1
        machines_A += 1
        return True
    return False

def method_5():
    global parts_X, parts_Y, parts_Z
    if parts_X >= 1 and parts_Y >= 1:
        parts_X -= 1
        parts_Y -= 1
        parts_Z += 1
        return True
    return False

# Cycle through the methods
methods = [method_1, method_2, method_3, method_4, method_5]
while any(method() for method in methods):
    pass

# Output the final quantities
result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
print(result)