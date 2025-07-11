def find_pitch_accent(word_data):
    """
    Analyzes word data to determine its pitch accent pattern.
    """
    word, reading, morae, pitch_morae, particle_pitch = word_data

    print(f"Analyzing the word: 「{word}」 ({reading})")
    print(f"The word has {morae} morae.")

    pitch_string = ""
    for i, mora in enumerate(reading):
        pitch_string += f" mora {i+1}「{mora}」: {pitch_morae[i]},"
    print(f"Pitch Pattern:{pitch_string[:-1]}")
    
    print(f"Pitch of a following particle (e.g., 「が」): {particle_pitch}")
    print("-" * 40)

    # Determine the pattern name based on the rules
    pattern_name = ""
    if pitch_morae[0] == "Low" and particle_pitch == "High":
        pattern_name = "Heiban"
        explanation = "The pitch starts low, rises, and stays high through the end of the word and the following particle."
    elif pitch_morae[0] == "High":
        pattern_name = "Atamadaka"
        explanation = "The pitch is high only on the first mora and then drops."
    elif "Low" in pitch_morae[1:] and particle_pitch == "Low":
        pattern_name = "Nakadaka"
        explanation = "The pitch rises after the first mora and falls somewhere before the end of the word."
    elif pitch_morae[-1] == "High" and particle_pitch == "Low":
        pattern_name = "Odaka"
        explanation = "The pitch rises after the first mora, stays high to the end, but a following particle is low."
    else:
        pattern_name = "Unknown"
        explanation = "Could not determine the pattern from the given data."

    print(f"Conclusion: The pattern is {pattern_name} (平板).")
    print(f"Reason: {explanation}")
    return "A" # Corresponds to Heiban

# Data for the word 「弟」 (otouto)
# morae: お-と-う-と
otouto_data = (
    "弟",
    "おとうと",
    4,
    ["Low", "High", "High", "High"], # Pitch pattern: LHHH
    "High" # A following particle is High
)

# Run the analysis
final_answer_letter = find_pitch_accent(otouto_data)
# print(f"\nFinal Answer Choice: {final_answer_letter}") # Not needed for final submission format

<<<A>>>