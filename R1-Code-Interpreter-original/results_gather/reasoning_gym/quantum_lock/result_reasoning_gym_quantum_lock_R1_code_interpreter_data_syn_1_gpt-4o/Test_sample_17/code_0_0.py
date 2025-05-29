def press_buttons(start, target):
    number = start
    light = 'red'
    sequence = []
    
    while number != target:
        if light == 'red':
            if number < target:
                number += 1
                sequence.append('C')
            else:
                number -= 1
                sequence.append('B')
        else:  # light is green
            if number * 2 <= target:
                number *= 2
                sequence.append('A')
            else:
                number += 1
                sequence.append('C')
        
        # Toggle light
        light = 'green' if light == 'red' else 'red'
    
    return ' → '.join(sequence)

print(press_buttons(0, 9))