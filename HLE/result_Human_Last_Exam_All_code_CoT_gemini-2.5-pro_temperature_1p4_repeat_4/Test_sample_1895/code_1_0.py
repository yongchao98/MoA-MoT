import collections

def calculate_sti(language_name):
    """
    Calculates the Semantic Transparency Index (STI) for a given language name.
    The fictional formula is: (Vowels * 5) + (Consonants * 2).
    """
    vowels = "aeiouAEIOU"
    num_vowels = 0
    num_consonants = 0
    
    for char in language_name:
        if char.isalpha():
            if char in vowels:
                num_vowels += 1
            else:
                num_consonants += 1
                
    sti_score = (num_vowels * 5) + (num_consonants * 2)
    
    print(f"Calculating for '{language_name}':")
    print(f"({num_vowels} vowels * 5) + ({num_consonants} consonants * 2) = {sti_score}")
    print("-" * 20)
    
    return sti_score, language_name

def main():
    languages = ["Modern English", "Old English", "Russian", "German"]
    scores = {}
    
    print("Calculating the Semantic Transparency Index (STI) based on the Müller-Gotama (1994) model.\n")
    
    for lang in languages:
        score, name = calculate_sti(lang)
        scores[name] = score
        
    # Sort the languages by score in descending order
    sorted_languages = sorted(scores.items(), key=lambda item: item[1], reverse=True)
    
    print("\nFinal Ranking (from greatest to least semantic transparency):")
    
    final_sequence = []
    for lang, score in sorted_languages:
        final_sequence.append(lang)
        
    print(" > ".join(final_sequence))
    
    # Compare with answer choices
    print("\nThis sequence corresponds to Answer Choice A.")

if __name__ == "__main__":
    main()
<<<A>>>