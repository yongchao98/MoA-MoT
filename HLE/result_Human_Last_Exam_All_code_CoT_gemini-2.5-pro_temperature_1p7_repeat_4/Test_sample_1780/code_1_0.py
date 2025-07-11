import itertools

def count_valid_arrangements():
    """
    Calculates the number of ways the letters L, N, S, W can be arranged
    based on a specific tip connection rule.
    """
    # Step 1: Define the properties of each letter's start (leftmost) and end (rightmost) tips.
    # 'T' stands for a tip at the Top, 'B' stands for a tip at the Bottom.
    letter_tips = {
        'L': {'start': 'T', 'end': 'B'},
        'N': {'start': 'B', 'end': 'T'},
        'S': {'start': 'T', 'end': 'B'},
        'W': {'start': 'T', 'end': 'T'}
    }
    
    letters = ['L', 'N', 'S', 'W']
    valid_arrangement_count = 0

    # Step 2: Define the connection rule function.
    # A connection is valid if the end tip of the first letter and the start tip
    # of the second letter are at different vertical levels (T-B or B-T).
    def can_connect(letter1, letter2):
        return letter_tips[letter1]['end'] != letter_tips[letter2]['start']

    # Step 3: Generate and check all permutations.
    all_permutations = list(itertools.permutations(letters))
    
    for p in all_permutations:
        # A permutation is a sequence like ('L', 'S', 'W', 'N').
        # Check if connections are valid for all adjacent letters in the sequence.
        if (can_connect(p[0], p[1]) and
            can_connect(p[1], p[2]) and
            can_connect(p[2], p[3])):
            
            # If all connections are valid, it's a valid arrangement.
            valid_arrangement_count += 1
            
    # Step 4: Output the final count.
    print(valid_arrangement_count)

count_valid_arrangements()
<<<2>>>