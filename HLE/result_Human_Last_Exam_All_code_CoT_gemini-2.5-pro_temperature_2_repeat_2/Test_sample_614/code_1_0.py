import itertools

# This script solves the riddle by programmatically checking all 120 possible
# arrangements of the numbers {0, 1, 2, 3, 4} against the rules.

# The key insight is that "strongest" and "weakest" refer to the number of
# letters in the English spelling of the numbers' names.
names = {0: "zero", 1: "one", 2: "two", 3: "three", 4: "four"}
strengths = {n: len(names[n]) for n in names}
min_strength = min(strengths.values())
max_strength = max(strengths.values())

# We check every possible permutation of the five numbers.
for p in itertools.permutations([0, 1, 2, 3, 4]):
    
    # Rule 1: "Number 1 is the first... and is the weakest"
    # The number 1 must be at the first position (index 0) and its strength
    # (number of letters) must be the minimum.
    if not (p[0] == 1 and strengths[1] == min_strength):
        continue

    # Rule 2: "...the fifth, who never had a thing and lastly follows"
    # The number 0 must be at the last position (index 4).
    if not (p[4] == 0):
        continue
        
    # Rule 3: "Number three is the strongest and it protects the last"
    # The number 3 must be the strongest and be in the second-to-last
    # position (index 3).
    if not (p[3] == 3 and strengths[3] == max_strength):
        continue

    # Rule 4: "Number 4... always goes behind [Number 2]"
    # The index of 4 must be one greater than the index of 2.
    if not p.index(4) == p.index(2) + 1:
        continue

    # If a permutation passes all checks, it's the solution.
    # The prompt asks to "output each number in the final equation!".
    # We will print each number of the solved sequence.
    print("The numbers in the final correct order are:")
    for number in p:
        print(number)
    
    # We can stop since there is only one unique solution.
    break