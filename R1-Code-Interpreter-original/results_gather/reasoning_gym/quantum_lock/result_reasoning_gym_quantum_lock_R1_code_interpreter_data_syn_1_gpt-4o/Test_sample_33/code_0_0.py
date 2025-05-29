def simulate_button_presses():
    start = 0
    target = 23
    light = 'red'
    sequence = []
    
    while start != target:
        if light == 'red':
            if start + 3 <= target:
                start += 3
                sequence.append('A')
            else:
                start -= 3
                sequence.append('C')
        else:  # light is green
            if start + 2 <= target:
                start += 2
                sequence.append('B')
            else:
                start += 3
                sequence.append('A')
        
        # Toggle light
        light = 'green' if light == 'red' else 'red'
    
    return ' → '.join(sequence)

print(simulate_button_presses())