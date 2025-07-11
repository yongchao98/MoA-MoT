# A list of the new and original sentences that are determined to be both true and meaningful.
true_and_meaningful_sentences = [
    "Satyendranath Bose baptized Bosons.",
    "Satyendranath Bose is a Bengalese.",
    "The conversation turned playful.",
    "Sushi contains fish.",
    "Our Tuna is playful.",
    "Our Tuna is a Bengalese."
]

# Calculate the total count of these sentences.
total_count = len(true_and_meaningful_sentences)

print("The true and meaningful sentences are:")
for sentence in true_and_meaningful_sentences:
    print(f"- {sentence}")

# Build the equation string dynamically to show each number.
equation_parts = ["1"] * total_count
equation_string = " + ".join(equation_parts)

print("\nFinal Calculation:")
# The final output is an equation where each '1' represents one true and meaningful sentence.
print(f"{equation_string} = {total_count}")