def explain_bayati_modulation():
    """
    This script explains the music theory behind the most common modulation
    in a Taqsim for Maqam Bayati on D and identifies the correct answer from a list.
    """
    # Maqam Bayati on D structure
    tonic = "D"
    second_degree = "E-half-flat"
    jins_bayati_on_d = ["D", second_degree, "F", "G"]

    # The most common modulation target
    correct_modulation_jins = "Saba"
    modulation_pivot_note = "E"  # Common name for the E-half-flat in this context
    correct_answer_choice = "H"

    print("Step-by-step analysis of the most common modulation from Maqam Bayati on D:")
    print("-" * 70)

    # Step 1: Define the base Maqam
    print("1. Understanding Maqam Bayati on D:")
    print(f"   - The tonic (qarar) is {tonic}.")
    print(f"   - The root scale fragment (Jins) is Jins Bayati, which consists of the notes: {', '.join(jins_bayati_on_d)}.")
    print(f"   - The most characteristic note is the second degree, '{second_degree}', which gives the maqam its unique sound.")
    print()

    # Step 2: Explain the most common modulation path
    print("2. Identifying the Common Modulation:")
    print("   - In a taqsim (improvisation), performers often modulate to create variety and emotional depth.")
    print("   - A classic and very common modulation in Bayati is to pivot on its unique second degree.")
    print(f"   - The performer treats this '{second_degree}' as a temporary tonic for a new Jins.")
    print()

    # Step 3: Evaluate the correct option
    print(f"3. Analyzing Answer Choice {correct_answer_choice}: Move to Jins {correct_modulation_jins} on {modulation_pivot_note}")
    print(f"   - This choice points directly to the classic modulation described above.")
    print(f"   - The 'E' in the answer refers to the {second_degree} of the Bayati scale.")
    print(f"   - When modulating to Jins {correct_modulation_jins}, the G-natural from Bayati is temporarily lowered to G-flat.")
    print("   - This creates the powerful and poignant sound of Saba, which contrasts beautifully with Bayati before resolving back.")
    print()

    # Step 4: Conclusion
    print("4. Conclusion:")
    print("   - While other modulations are theoretically possible, the move to Jins Saba on the second degree is a hallmark of Bayati performance practice.")
    print("   - The other options listed are either very rare or structurally distant, making them highly unusual choices for a typical Bayati taqsim.")
    print("-" * 70)
    print("Therefore, the most common modulation listed is to Jins Saba on E.")

explain_bayati_modulation()
print("<<<H>>>")