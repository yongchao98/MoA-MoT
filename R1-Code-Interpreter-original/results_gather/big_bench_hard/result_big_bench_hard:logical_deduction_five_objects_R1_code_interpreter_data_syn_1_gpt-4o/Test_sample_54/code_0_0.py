from itertools import permutations

books = ['orange', 'gray', 'purple', 'yellow', 'green']
positions = list(permutations(books))

for pos in positions:
    if pos[4] == 'green':  # Green book is the rightmost
        if pos.index('gray') > pos.index('orange'):  # Gray is to the right of Orange
            if pos.index('purple') < pos.index('yellow'):  # Purple is to the left of Yellow
                if pos.index('purple') > pos.index('gray'):  # Purple is to the right of Gray
                    print(pos)
                    break