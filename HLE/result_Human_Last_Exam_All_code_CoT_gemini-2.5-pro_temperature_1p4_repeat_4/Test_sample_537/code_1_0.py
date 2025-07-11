import collections

def num_to_words(n):
    """Converts a number to its English word representation."""
    # A simple dictionary for the numbers needed, to avoid external libraries.
    d = {
        0: 'zero', 1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five',
        6: 'six', 7: 'seven', 8: 'eight', 9: 'nine', 10: 'ten',
        11: 'eleven', 12: 'twelve', 13: 'thirteen', 14: 'fourteen', 15: 'fifteen',
        16: 'sixteen', 17: 'seventeen', 18: 'eighteen', 19: 'nineteen', 20: 'twenty',
        21: 'twentyone', 22: 'twentytwo', 23: 'twentythree', 24: 'twentyfour',
        25: 'twentyfive', 26: 'twentysix', 27: 'twentyseven', 28: 'twentyeight',
        29: 'twentynine', 30: 'thirty', 31: 'thirtyone', 32: 'thirtytwo',
        33: 'thirtythree', 34: 'thirtyfour', 35: 'thirtyfive', 36: 'thirtysix',
        37: 'thirtyseven', 38: 'thirtyeight', 39: 'thirtynine', 40: 'forty'
    }
    # For this puzzle, we only need to represent numbers as a single word.
    # We remove spaces and hyphens if they were present.
    return d[n].replace(' ', '').replace('-', '') if n in d else ''

def calculate_match_score(word, number_as_word):
    """Calculates a normalized score for how well a word matches a number's name."""
    if not number_as_word:
        return 0
    word_counts = collections.Counter(word)
    number_word_counts = collections.Counter(number_as_word)
    
    # Calculate the number of common letters
    common_letters_count = sum((word_counts & number_word_counts).values())
    
    # Normalize the score by the length of the number's name
    return common_letters_count / len(number_as_word)

def find_best_sequence():
    """Searches for the number sequence that best matches the given words."""
    words = ['utpshtheas', 'fkuhu', 'numhkatnatu', 'anhaeketn']
    
    best_score = -1
    best_sequence = None
    
    # Iterate through possible starting numbers for the Fibonacci-like sequence
    for n1 in range(1, 26):
        for n2 in range(1, 26):
            if n1 == n2:
                continue
            
            n3 = n1 + n2
            n4 = n2 + n3
            
            sequence = [n1, n2, n3, n4]
            number_words = [num_to_words(n) for n in sequence]
            
            if any(not s for s in number_words):
                continue
                
            current_total_score = 0
            for i in range(4):
                current_total_score += calculate_match_score(words[i], number_words[i])
            
            if current_total_score > best_score:
                best_score = current_total_score
                best_sequence = sequence
                
    return best_sequence

# Find the sequence and print the final result
final_sequence = find_best_sequence()

if final_sequence:
    n2, n3, n4 = final_sequence[1], final_sequence[2], final_sequence[3]
    # The last word 'anhaeketn' decodes to the last number in the sequence.
    # The puzzle asks for the final equation.
    print(f"The decoded last word is {n4}.")
    print("The final equation is:")
    print(f"{n2} + {n3} = {n4}")
else:
    print("Could not find a matching sequence.")
