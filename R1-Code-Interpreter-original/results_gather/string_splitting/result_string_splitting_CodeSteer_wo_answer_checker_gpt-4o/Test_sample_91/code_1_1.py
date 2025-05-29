# Initialize the counts of machines and parts in a dictionary
state = {
    'machines_A': 4,
    'machines_B': 4,
    'machines_C': 6,
    'parts_X': 0,
    'parts_Y': 0,
    'parts_Z': 0
}

# Define the methods as functions
def method_1(state):
    if state['machines_A'] >= 1:
        state['machines_A'] -= 1
        state['parts_X'] += 2
        state['parts_Y'] += 1
        return True
    return False

def method_2(state):
    if state['machines_B'] >= 2:
        state['machines_B'] -= 2
        state['parts_X'] += 1
        return True
    return False

def method_3(state):
    if state['machines_C'] >= 2:
        state['machines_C'] -= 2
        state['parts_Y'] += 1
        return True
    return False

def method_4(state):
    if state['machines_B'] >= 1 and state['machines_C'] >= 1:
        state['machines_B'] -= 1
        state['machines_C'] -= 1
        state['machines_A'] += 1
        return True
    return False

def method_5(state):
    if state['parts_X'] >= 1 and state['parts_Y'] >= 1:
        state['parts_X'] -= 1
        state['parts_Y'] -= 1
        state['parts_Z'] += 1
        return True
    return False

# List of methods to apply in order
methods = [method_1, method_2, method_3, method_4, method_5]

# Function to execute the methods in sequence
def execute_methods(state):
    while True:
        executed = False
        for method in methods:
            if method(state):
                executed = True
                break
        if not executed:
            break

# Execute the methods
execute_methods(state)

# Output the final counts
result = [str(state['machines_A']), str(state['machines_B']), str(state['machines_C']),
          str(state['parts_X']), str(state['parts_Y']), str(state['parts_Z'])]
print(result)