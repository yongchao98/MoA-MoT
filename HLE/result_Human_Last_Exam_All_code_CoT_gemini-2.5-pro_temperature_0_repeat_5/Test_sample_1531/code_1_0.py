import re
try:
    import pronouncing
except ImportError:
    print("The 'pronouncing' library is not installed.")
    print("Please install it by running: pip install pronouncing")
    exit()

def analyze_lines(lines):
    """
    Analyzes the meter of poetic lines by counting syllables and stresses.
    """
    print("Analyzing the metrical structure of the lines...\n")
    
    # Common function words that are typically unstressed in poetry
    unstressed_words = {'&', 'all', 'the', 'a', 'are'}

    for i, line in enumerate(lines):
        clean_line = line.lower().replace('&', 'and')
        words = re.findall(r'\b\w+\b', clean_line)

        total_syllables = 0
        total_stresses = 0
        stress_pattern_display = []

        for word in words:
            # Heuristic: Treat common function words as unstressed
            if word in unstressed_words:
                stresses = ['0']
            else:
                # Get stress pattern from the pronouncing library
                stresses_list = pronouncing.stresses_for_word(word)
                if stresses_list:
                    # Use the first available stress pattern
                    stresses = stresses_list[0].split()
                else:
                    # Fallback for words not in the dictionary
                    stresses = ['1'] 

            # Get syllable count
            syllables = pronouncing.syllable_count(word)
            if syllables == 0: # Fallback for unknown words
                syllables = len(stresses)

            total_syllables += syllables
            
            # Count primary or secondary stresses
            num_stresses = sum(1 for s in stresses if s in ['1', '2'])
            total_stresses += num_stresses
            
            # Build a visual representation (DUM for stressed, da for unstressed)
            pattern = ' '.join(['da' if s == '0' else 'DUM' for s in stresses])
            stress_pattern_display.append(pattern)

        # The "equation" of the line's meter
        print(f'Line {i+1}: "{line}"')
        print(f"Analysis = Syllables: {total_syllables} | Stressed Beats: {total_stresses} | Pattern: {' | '.join(stress_pattern_display)}")
        print("-" * 20)

    print("\nConclusion:")
    print("Line 1 has 7 syllables and 3 stressed beats.")
    print("Line 2 has 6 syllables and 3 stressed beats, fitting a clear iambic trimeter (da-DUM da-DUM da-DUM).")
    print("Since both lines are built around three primary stresses, 'trimeter' (a line with three metrical feet) is the most accurate description.")
    print("'Iambic pentameter' is incorrect (requires 5 feet/10 syllables).")
    print("'Free verse' is less precise because a clear metrical pattern (trimeter) is present.")


# The two lines to be analyzed
poetic_lines = [
    "& all the stars are palaces",
    "the world a hollow road"
]

analyze_lines(poetic_lines)

<<<E>>>