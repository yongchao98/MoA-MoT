# Possible numbers and letters based on the analysis
numbers = [2, 8]
letters = ['M', 'P']

# Check all possible combinations of the letters
for z in letters:
    for w in letters:
        if z != w:
            # Check the conditions
            if (z == 'M' and w == 'P') or (z == 'P' and w == 'M'):
                print([numbers[0], numbers[1], z, w])