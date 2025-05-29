def simulate_button_presses():
    number = 0
    light = 'red'
    sequence = []

    while number != 9:
        if light == 'red':
            # Press B
            number += 1
            light = 'green'
            sequence.append('B')
        elif light == 'green':
            # Press A
            number += 2
            light = 'red'
            sequence.append('A')

    return ' → '.join(sequence)

print(simulate_button_presses())