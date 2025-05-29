from itertools import permutations

golfers = ['Joe', 'Dan', 'Ada', 'Amy', 'Rob', 'Mya', 'Mel']
constraints = [
    lambda order: order.index('Amy') == 2,  # Amy finished third
    lambda order: order.index('Mel') == 4,  # Mel finished third-to-last
    lambda order: order.index('Ada') > order.index('Amy'),  # Ada below Amy
    lambda order: order.index('Joe') > order.index('Dan'),  # Joe below Dan
    lambda order: order.index('Dan') > order.index('Ada'),  # Dan below Ada
    lambda order: order.index('Rob') > order.index('Mya')  # Rob below Mya
]

for order in permutations(golfers):
    if all(constraint(order) for constraint in constraints):
        print(order)
        break