# Possible numbers and letters
numbers = [3, 4, 6, 7]
letters = ['A', 'B', 'C', 'D', 'E']

# Check all combinations
for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:
                        # Check each condition
                        if (n1 == 4 and n2 == 6 and l1 == 'T' and l2 == 'E'):
                            print([n1, n2, l1, l2])