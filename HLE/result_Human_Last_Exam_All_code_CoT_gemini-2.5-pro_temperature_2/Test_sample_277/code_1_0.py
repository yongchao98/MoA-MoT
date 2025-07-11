from fractions import Fraction

# Step 1: Define initial parameters.
# The song starts with note 'G', and its frequency is given as 392 Hz.
initial_frequency = 392

# The song ends with the note 'C'.
start_note_name = 'G'
end_note_name = 'C'

# Step 2: Define the just intonation ratios for the start and end notes.
# These ratios are defined relative to a tonic. Since the notes are in C-major,
# we use C as the tonic with a ratio of 1/1.
just_ratios = {
    'C': Fraction(1, 1),  # Tonic
    'G': Fraction(3, 2)   # Perfect Fifth
}

# Step 3: Calculate the cumulative frequency ratio.
# Because the frequency of each note is calculated relative to the one just
# before it, the cumulative effect is a telescoping product. The overall ratio
# is the just ratio of the end note divided by the just ratio of the start note.
total_ratio = just_ratios[end_note_name] / just_ratios[start_note_name]

# Step 4: Calculate the final frequency.
# We multiply the initial frequency by the total cumulative ratio.
final_frequency = Fraction(initial_frequency) * total_ratio

# Step 5: Extract the numerator (a) and denominator (b).
a = final_frequency.numerator
b = final_frequency.denominator

# Step 6: Print the final output as requested.
# The instruction asks to output each number in the final equation.
initial_freq_val = initial_frequency
total_ratio_num = total_ratio.numerator
total_ratio_den = total_ratio.denominator

print(f"Starting with G at {initial_freq_val} Hz and ending at C, the final frequency calculation is:")
print(f"{initial_freq_val} * ({total_ratio_num}/{total_ratio_den}) = {a}/{b} Hertz")
print("The final answer [a, b], where a and b are coprime, is:")
print([a, b])