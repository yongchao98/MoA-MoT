def simulate_sequence():
    # Initial state
    number = 0
    light = 'red'
    sequence = []

    # Define button operations
    def press_A(num, light):
        return num * 2, 'green' if light == 'red' else 'red'

    def press_B(num, light):
        if light == 'red':
            return num * 2, 'green'
        return num, light

    def press_C(num, light):
        if light == 'green':
            return num + 1, 'red'
        return num, light

    # Simulate the sequence
    number, light = press_C(number, light)
    sequence.append('C')
    number, light = press_A(number, light)
    sequence.append('A')
    number, light = press_A(number, light)
    sequence.append('A')
    number, light = press_A(number, light)
    sequence.append('A')

    return sequence

# Execute the simulation
sequence = simulate_sequence()
print(" → ".join(sequence))