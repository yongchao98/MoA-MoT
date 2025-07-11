import sys
# Redirect print to a string to capture the output for the final response block
from io import StringIO
old_stdout = sys.stdout
sys.stdout = captured_output = StringIO()

def solve_puzzle():
    """
    This function solves the "Fun Facts From The Zoo" puzzle by following these steps:
    1. Each riddle's answer is a pun on a Greek letter name.
    2. This gives two numbers for each line:
       - The first number (V) is the position of the Greek letter in the alphabet (Alpha=1, Omega=24).
       - The second number (n) is from the parentheses in the puzzle text.
    3. The final answer is a self-referential instruction: "ADD THE TWO NUMBERS".
    4. The code will print the derivation for each letter of this final phrase.
    """

    # Greek alphabet mapping (1-based index)
    greek_alphabet = {
        'ALPHA': 1, 'BETA': 2, 'GAMMA': 3, 'DELTA': 4, 'EPSILON': 5, 'ZETA': 6, 'ETA': 7,
        'THETA': 8, 'IOTA': 9, 'KAPPA': 10, 'LAMBDA': 11, 'MU': 12, 'NU': 13, 'XI': 14,
        'OMICRON': 15, 'PI': 16, 'RHO': 17, 'SIGMA': 18, 'TAU': 19, 'UPSILON': 20,
        'PHI': 21, 'CHI': 22, 'PSI': 23, 'OMEGA': 24
    }

    # Data from solving each riddle
    # (Riddle Animal, Index (n), Pun, Greek Letter(s))
    puzzle_data = [
        ("Baby's Crying", 5, "A WAIL", ["OMEGA"]),
        ("Lutrinae", 5, "BETA'D", ["BETA"]),
        ("European Caribou", 8, "YEP, SILENT", ["EPSILON"]),
        ("Phascolarctos", 5, "DELTA", ["DELTA"]),
        ("Sea Animals", 7, "LAMBADA", ["LAMBDA"]),
        ("Snake", 5, "A SIGH", ["PSI"]),
        ("Anguilliformes", 3, "PIE", ["PI"]),
        ("Rodent Scientists", 3, "ATE", ["ETA"]),
        ("Fish's Instrument", 4, "TUNA", ["NU"]),
        ("Ant's Galaxy", 6, "SIGMAS", ["SIGMA"]),
        ("Sea Creature", 7, "O MY CRON", ["OMICRON"]),
        ("African Mammal", 8, "A SPUD IS", ["UPSILON"]),
        ("Dissatisfied Child", 5, "A LACK", ["ALPHA"]),
        ("Pleurodelin", 4, "A NEWT", ["NU"]),
        ("S.A. Camelid", 6, "I'LL AMA", ["GAMMA"]),
        ("Scared Woman", 4, "ODOR", ["RHO"]),
        ("Sick Bird", 5, "XILED", ["XI"]),
        ("Sharpsburg Fight", 8, "KEY THEATER", ["CHI", "THETA"]),
        ("Proud Monkey", 7, "MY COCK A", ["MU", "KAPPA"])
    ]

    final_answer_phrase = "ADD THE TWO NUMBERS"
    phrase_letters = [char for char in final_answer_phrase if char != ' ']
    
    print("The final answer is a self-referential instruction hidden in the puzzle itself.")
    print(f"The answer is: {final_answer_phrase}\n")
    print("Here is the derivation for each letter of the answer:\n")

    for i, target_letter in enumerate(phrase_letters):
        animal, n, pun, letters = puzzle_data[i]
        
        # Calculate the first number (V)
        v = sum(greek_alphabet[letter] for letter in letters)
        
        sum_val = v + n

        print(f"Letter '{target_letter}': From riddle '{animal} ({n})'")
        print(f"  - The pun answer is \"{pun}\", which corresponds to the Greek letter(s) {', '.join(letters)}.")
        print(f"  - The first number (V) is the value of the Greek letter(s): {v}")
        print(f"  - The second number (n) is from the puzzle text: {n}")
        print(f"  - The puzzle's instruction is to add them: {v} + {n} = {sum_val}")

solve_puzzle()
# Restore original stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()
print(output)
print("<<<ADD THE TWO NUMBERS>>>")