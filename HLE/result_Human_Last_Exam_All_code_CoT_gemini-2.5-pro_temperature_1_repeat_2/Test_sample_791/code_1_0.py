def generate_adjective_sequence(n):
    """
    Generates a sequence of 'red' adjectives.
    The grammar's left recursion REC = REC ADJ allows for such sequences.
    """
    if n <= 0:
        return ""
    # Using 'red' as a stand-in for the repeating ADJ
    return " ".join(["red"] * n)

def generate_sentence(num_adjectives):
    """
    Generates a valid sentence from the language L using the rule S -> N ADJ.
    This demonstrates that sentences can be arbitrarily long.
    """
    # S -> N ADJ
    noun = "snakes"
    adjectives = generate_adjective_sequence(num_adjectives)
    subject = f"{noun} {adjectives}"

    # V -> 'swim'
    verb = "swim"
    # C -> 'well'
    complement = "well"

    # L -> S V C '.'
    return f"{subject} {verb} {complement}."

print("Demonstrating that the language L has no longest sentence:")
print("A valid sentence in L is: 'red frogs swim swiftly.'")
print("-" * 20)
print("Here are other valid sentences of increasing length:")
for i in range(1, 5):
    print(f"Sentence with {i} adjective(s): {generate_sentence(i)}")
