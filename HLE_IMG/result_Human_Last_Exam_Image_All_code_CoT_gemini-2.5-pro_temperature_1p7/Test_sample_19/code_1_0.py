class MusicTheoryAnalyzer:
    """
    A class to analyze the provided musical excerpt and answer the user's questions.
    """
    def __init__(self):
        self.home_key = "C# minor"
        self.relative_major = "E Major"

    def explain_analysis(self):
        """
        Prints the step-by-step analysis for the three-part question.
        """
        print("Here is a step-by-step analysis of the excerpt from Beethoven's 'Moonlight' Sonata:")
        print("="*80)

        # Part 1: Identify the key of modulation
        print("\nPart 1: To what key the music modulates in measure 11 and the beginning of measure 12?")
        print("-" * 80)
        print(f"1. The piece's home key, based on the four-sharp key signature and opening harmony, is {self.home_key}.")
        print("2. In measure 11, the left hand arpeggiates the notes B, D#, and F#. These three notes form a B Major chord.")
        print("3. On the first beat of measure 12, both hands play the notes of a B Major chord (B-D#-F#).")
        print("\nAnswer 1: The music modulates to the key of B Major.")
        print("="*80)
        
        # Part 2: Justify the modulation
        print("\nPart 2: What is the justification for this particular modulation?")
        print("-" * 80)
        print(f"1. The home key is {self.home_key}. Its relative major is {self.relative_major}.")
        print(f"2. The preceding measures (9-10) clearly establish the key of {self.relative_major}, which is the mediant (III) chord of {self.home_key}.")
        print(f"3. B Major is the dominant (the V chord) of {self.relative_major}.")
        print("\nAnswer 2: The justification is that the music has modulated to the dominant of the relative major. This is a common and strong harmonic progression used to extend the section in the relative major key.")
        print("="*80)

        # Part 3: Provide the Roman Numeral
        print("\nPart 3: What is the complete Roman numeral for the first beat of measure 11?")
        print("-" * 80)
        print("1. Roman numeral analysis relates chords back to the original home key.")
        print(f"2. The chord is B Major, and the home key is {self.home_key}.")
        print("3. As established, B Major is functioning as the Dominant of the Mediant chord (E Major).")
        print("4. This 'Dominant of the Mediant' function is written as a secondary dominant.")
        print("5. The Mediant is Roman numeral III (for the 3rd degree of the scale).")
        print("6. The Dominant is Roman numeral V (for the 5th degree of the scale).")
        print("\nAnswer 3: The complete and accurate Roman numeral marking is V/III.")
        print("="*80)

# Create an instance and run the analysis
analyzer = MusicTheoryAnalyzer()
analyzer.explain_analysis()