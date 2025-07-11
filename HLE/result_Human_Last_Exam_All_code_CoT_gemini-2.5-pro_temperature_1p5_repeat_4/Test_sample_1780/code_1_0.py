import itertools

def solve_letter_arrangement():
    """
    Calculates the number of ways to arrange the letters "L", "N", "S", "W"
    based on a specific connection rule.
    """
    # Step 1: Define the properties of each letter's unconnected tips.
    # T = Top, B = Bottom, E = End of a horizontal line.
    letter_tips = {
        'L': {'left': 'T', 'right': 'E'},
        'N': {'left': 'T', 'right': 'B'},
        'S': {'left': 'T', 'right': 'B'},
        'W': {'left': 'T', 'right': 'T'}
    }

    # Step 2: Define the connection logic based on the prompt's examples.
    # "S" (right: B) can connect to "W" (left: T) => B -> T is valid.
    # "W" (right: T) cannot connect to "S" (left: T) => T -> T is invalid.
    # We infer that L's end tip can connect to a top tip (E -> T is valid).
    def can_connect(letter1, letter2):
        """Checks if the right tip of letter1 can connect to the left of letter2."""
        right_tip1 = letter_tips[letter1]['right']
        left_tip2 = letter_tips[letter2]['left']

        if right_tip1 == 'B' and left_tip2 == 'T':
            return True
        if right_tip1 == 'E' and left_tip2 == 'T':
            return True
        # Any other connection, including T -> T, is invalid.
        return False

    # Step 3: Generate all permutations of the letters.
    letters = ['L', 'N', 'S', 'W']
    all_permutations = list(itertools.permutations(letters))
    
    valid_arrangement_count = 0
    valid_arrangements = []

    # Step 4: Check each permutation for valid sequential connections.
    for p in all_permutations:
        # A valid arrangement requires all adjacent letters to connect properly.
        # p[0]->p[1], p[1]->p[2], and p[2]->p[3]
        if (can_connect(p[0], p[1]) and
            can_connect(p[1], p[2]) and
            can_connect(p[2], p[3])):
            valid_arrangement_count += 1
            valid_arrangements.append("".join(p))
            
    # As an interesting note, we can deduce logically that 'W' must always
    # be the last letter, as it cannot connect to any other letter.
    # The problem then becomes finding the number of permutations of ('L', 'N', 'S'),
    # which is 3! = 6. Let's confirm this with our code's findings.

    # print("The valid arrangements are:")
    # for arr in valid_arrangements:
    #     print(arr)
    # print("\nTotal number of ways:")
    print(valid_arrangement_count)


solve_letter_arrangement()