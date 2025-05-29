from itertools import permutations

birds = ['cardinal', 'hummingbird', 'raven', 'robin', 'crow', 'quail', 'blue jay']
positions = list(permutations(birds))

for pos in positions:
    if (pos[0] == 'cardinal' and
        pos[3] == 'hummingbird' and
        pos[4] == 'blue jay' and
        pos.index('hummingbird') < pos.index('crow') and
        pos.index('crow') < pos.index('quail') and
        pos.index('raven') < pos.index('robin')):
        print(pos)
        break