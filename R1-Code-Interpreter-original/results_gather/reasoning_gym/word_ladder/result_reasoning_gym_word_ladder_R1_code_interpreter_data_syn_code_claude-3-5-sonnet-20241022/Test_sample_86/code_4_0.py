# Let's add debug printing to verify our word list and transitions
words = {'LAST', 'LOST', 'COST', 'COSH', 'COHO'}

def is_one_letter_diff(word1, word2):
    diff = 0
    for c1, c2 in zip(word1, word2):
        if c1 != c2:
            diff += 1
    return diff == 1

print("Word list:", words)
print("\nVerifying transitions:")
transitions = [('LAST', 'LOST'), ('LOST', 'COST'), ('COST', 'COSH'), ('COSH', 'COHO')]
for w1, w2 in transitions:
    print(f"{w1} -> {w2}: {'Valid' if is_one_letter_diff(w1, w2) and w1 in words and w2 in words else 'Invalid'}")
    if not is_one_letter_diff(w1, w2):
        print(f"  Not one letter difference")
    if w1 not in words:
        print(f"  {w1} not in word list")
    if w2 not in words:
        print(f"  {w2} not in word list")