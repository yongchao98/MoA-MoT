def reach_target(start, target):
    sequence = []
    number = start
    light = 'red'  # initial light state

    while number != target:
        if light == 'red':
            if number < target:
                sequence.append('A')
                number += 1
                light = 'green'
            else:
                sequence.append('C')
                number -= 1
                light = 'green'
        else:  # light is green
            if number * 2 <= target:
                sequence.append('B')
                number *= 2
                light = 'red'
            else:
                sequence.append('A')
                number += 1
                light = 'red'

    return ' → '.join(sequence)

print(reach_target(0, 9))