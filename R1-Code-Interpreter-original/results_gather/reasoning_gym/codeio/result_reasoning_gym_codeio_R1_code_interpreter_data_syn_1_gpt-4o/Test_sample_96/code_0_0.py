from collections import namedtuple

FIELDS = [
    "children",
    "cats",
    "samoyeds",
    "pomeranians",
    "akitas",
    "vizslas",
    "goldfish",
    "trees",
    "cars",
    "perfumes"
]

SUE_DATA = dict(
    children=lambda x: x == 3,
    cats=lambda x: x > 7,
    samoyeds=lambda x: x == 2,
    pomeranians=lambda x: x < 3,
    akitas=lambda x: x == 0,
    vizslas=lambda x: x == 0,
    goldfish=lambda x: x < 5,
    trees=lambda x: x > 3,
    cars=lambda x: x == 2,
    perfumes=lambda x: x == 1
)

Sue = namedtuple("Sue", ["number"] + FIELDS, defaults=[None] * (len(FIELDS) + 1))

def calculate_score(sue):
    matching_fields = 0

    for i, field in enumerate(FIELDS, 1):
        if sue[i] is None:
            continue
        if SUE_DATA[field](sue[i]):
            matching_fields += 1

    return matching_fields

def parse(raw):
    sues = []

    for sue in raw:
        parts = sue.replace(":", "").replace(",", "").split()
        deets = {}
        i = 2
        while i < len(parts):
            deets[parts[i]] = int(parts[i + 1])
            i += 2

        sues.append(Sue(number=int(parts[1]), **deets))

    return sues

def main_solution(raw_data):
    sues = parse(raw_data)
    best_match = max(sues, key=lambda sue: calculate_score(sue))
    return best_match.number

# Test the function with the provided input
raw_data = [
    'Sue 1, vizslas: 0, samoyeds: 0, trees: 4, children: 2, cats: 1, cars: 3, akitas: 1, pomeranians: 3, perfumes: 2, goldfish: 3',
    'Sue 2, trees: 3, vizslas: 0, cars: 1, pomeranians: 4, perfumes: 1, samoyeds: 0, akitas: 0',
    'Sue 3, akitas: 3, trees: 2, goldfish: 3, samoyeds: 0, vizslas: 3, perfumes: 3, pomeranians: 4, cats: 4, children: 2',
    'Sue 4, goldfish: 3, cars: 0, perfumes: 2, cats: 0, trees: 10, pomeranians: 2, vizslas: 0, samoyeds: 1, children: 1',
    'Sue 5, cars: 1, cats: 10, perfumes: 2',
    'Sue 6, vizslas: 1, pomeranians: 0, cats: 10, goldfish: 3, children: 3, akitas: 2, trees: 0, cars: 0',
    'Sue 7, cats: 10, cars: 3, children: 3, akitas: 3, trees: 1, perfumes: 3, pomeranians: 3, goldfish: 5, vizslas: 2',
    'Sue 8, pomeranians: 0, cars: 0, perfumes: 1, akitas: 2',
    'Sue 9, pomeranians: 4, trees: 3, goldfish: 2, cats: 2, perfumes: 2, cars: 0, children: 3, vizslas: 2, akitas: 1'
]

print(main_solution(raw_data))