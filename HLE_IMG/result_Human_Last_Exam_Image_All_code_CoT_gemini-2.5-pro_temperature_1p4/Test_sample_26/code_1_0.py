import sys
# Redirect stdout to a variable to format the final output as requested.
original_stdout = sys.stdout
from io import StringIO
sys.stdout = captured_output = StringIO()


def analyze_chord():
    """
    Analyzes the chord in measure 30 of Beethoven's Pathetique Sonata, 1st mvt.
    and determines the correct Roman numeral.
    """
    key = "C minor"
    measure_number = 30
    notes_in_chord = ["G", "B-natural", "D", "F-sharp"]
    root_of_chord = "G"
    scale_degree = 5
    function = "Dominant"
    roman_numeral = "V"

    print(f"Musical Analysis for Measure {measure_number}:")
    print(f"1. The key of the piece is {key}.")
    print(f"2. The notes in measure {measure_number} are: {', '.join(notes_in_chord)}.")
    print(f"3. The root of the chord is {root_of_chord}, which is the fifth (dominant) scale degree of {key}.")
    print(f"4. The chord's function is {function}, as it leads strongly to the tonic chord (C minor) in the next measure.")
    print("5. The B-natural is the standard raised leading-tone for the key of C minor.")
    print("6. The F-sharp is a chromatic embellishment that strengthens the dominant pull.")
    print("\nBased on this analysis, the correct Roman numeral for the chord is:")
    print(roman_numeral)

analyze_chord()

# Get the content from captured_output and add the final answer tag.
final_output = captured_output.getvalue().strip()
sys.stdout = original_stdout
print(final_output)
final_answer = final_output.splitlines()[-1]
print(f'<<<{final_answer}>>>')
