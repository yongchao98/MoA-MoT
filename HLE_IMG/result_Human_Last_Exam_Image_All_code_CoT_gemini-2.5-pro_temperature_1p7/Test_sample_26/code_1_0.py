import sys

def analyze_beethoven_chord():
    """
    Provides a music theory analysis of the chord in measure 30
    of Beethoven's "Pathetique" Sonata, 1st movement, and prints
    the correct Roman numeral.
    """

    if sys.stdout.encoding != 'utf-8':
        sys.stdout.reconfigure(encoding='utf-8')

    print("Step 1: Identifying the musical context.")
    print("The piece is Beethoven's Piano Sonata No. 8, 'Pathetique'.")
    print("The first movement is in the key of C minor.")
    print("We are analyzing the chord in measure 30.\n")

    print("Step 2: Identifying the notes and chords.")
    print("The chord in measure 30 is composed of the notes C-sharp, E-natural, and G-sharp. This forms a C-sharp minor triad.")
    print("This chord resolves to the D minor triad (D-F-A) in the next measure (m. 31).\n")

    print("Step 3: Determining the function of the chords.")
    print("In the home key of C minor, the D minor chord is built on the 2nd scale degree. Therefore, it is the 'ii' chord.")
    print("The C-sharp minor chord in measure 30 directly precedes and resolves to the 'ii' chord, giving it a secondary function.\n")

    print("Step 4: Constructing the Roman Numeral.")
    print("The root of the chord, C-sharp, is the leading tone to D (the root of the 'ii' chord).")
    print("A chord built on a leading tone is labeled 'vii'. The diminished symbol '°' is standard for leading-tone chords.")
    print("Since this chord is applied to 'ii', it is a 'vii of ii', written with a slash.")
    print("Beethoven's specific spelling as a minor triad is an alteration for voice-leading, but the chord's function remains that of a secondary leading-tone chord.\n")

    print("Step 5: Final Answer Construction.")
    part1 = "vii"
    part2 = "\u00b0"
    part3 = "/"
    part4 = "ii"
    
    print("The Roman numeral is composed of the following parts:")
    print(f"- The leading-tone number: {part1}")
    print(f"- The diminished quality symbol: {part2}")
    print(f"- The 'applied to' separator: {part3}")
    print(f"- The target chord number: {part4}")

    final_answer = part1 + part2 + part3 + part4
    print(f"\nThus, the correct Roman numeral for measure 30 is: {final_answer}")

if __name__ == '__main__':
    analyze_beethoven_chord()
<<<vii°/ii>>>