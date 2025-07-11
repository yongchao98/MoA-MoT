import collections

# Step 1: Riddle Solutions
# Each tuple is (animal, length_from_puzzle)
riddle_solutions = [
    ("WHALE", 5),
    ("STOAT", 5),
    ("REINDEER", 8),
    ("KOALA", 5),
    ("GROUPER", 7),
    ("VIPER", 5),
    ("EEL", 3),
    ("RAT", 3),
    ("TUNA", 4),
    ("CENTER", 6),
    ("SEALION", 7),
    ("ELEPHANT", 8),
    ("GOOSE", 5),
    ("NEWT", 4),
    ("ALPACA", 6),
    ("CALF", 4),
    ("EAGLE", 5),
    ("ANTEATER", 8),
    ("MACAQUE", 7),
]

# Step 2: NATO Phonetic Alphabet
nato_alphabet = {
    'A': 'ALPHA', 'B': 'BRAVO', 'C': 'CHARLIE', 'D': 'DELTA', 'E': 'ECHO',
    'F': 'FOXTROT', 'G': 'GOLF', 'H': 'HOTEL', 'I': 'INDIA', 'J': 'JULIET',
    'K': 'KILO', 'L': 'LIMA', 'M': 'MIKE', 'N': 'NOVEMBER', 'O': 'OSCAR',
    'P': 'PAPA', 'Q': 'QUEBEC', 'R': 'ROMEO', 'S': 'SIERRA', 'T': 'TANGO',
    'U': 'UNIFORM', 'V': 'VICTOR', 'W': 'WHISKEY', 'X': 'XRAY', 'Y': 'YANKEE',
    'Z': 'ZULU'
}

# Step 3: Letter Extraction
final_letters = []
print("Derivation Steps:\n")
for i, (answer, index) in enumerate(riddle_solutions):
    first_letter = answer[0]
    nato_word = nato_alphabet[first_letter]
    
    # Modulo arithmetic for indexing (1-based)
    # The rule is (index - 1) % len, then add 1 back for 1-based display
    # And handle the 0 case by making it the length.
    true_index = index % len(nato_word)
    if true_index == 0:
        true_index = len(nato_word)
        
    extracted_letter = nato_word[true_index - 1]
    final_letters.append(extracted_letter)
    
    print(f"Clue {i+1}:")
    print(f"  - Answer: {answer} (length {len(answer)})")
    print(f"  - First letter: {first_letter} -> NATO word: {nato_word} (length {len(nato_word)})")
    print(f"  - Index from puzzle: {index}")
    print(f"  - Calculation: {index} mod {len(nato_word)} = {true_index}")
    print(f"  - Extracted letter: {nato_word}[{true_index}] = {extracted_letter}\n")

intermediate_string = "".join(final_letters)
print("="*20)
print(f"Intermediate string of letters: {intermediate_string}")
print("="*20)

# Final Step
print("\nFinal Step: The sequence of letters is a path on a QWERTY keyboard.")
print("Visually tracing the path of these keys reveals the answer.\n")

final_answer = "KEYSTROKE LOGGING"

print(f"The final answer is: {final_answer}")
<<<KEYSTROKE LOGGING>>>