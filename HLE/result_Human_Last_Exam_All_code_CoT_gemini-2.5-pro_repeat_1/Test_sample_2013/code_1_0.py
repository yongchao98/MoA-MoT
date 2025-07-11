# 1. Define the details of the musical piece.
# The joke described occurs in Camille Saint-Saëns's second piano concerto.
composer_surname = "Saint-Saens"
opus_number = 22

# 2. Define the measures for the musical joke.
# In the finale (Presto), the cello/bass flourish begins in measure 105.
# They are off-beat for measures 106 and 107.
# They skip a beat and recover in measure 108.
start_measure = 105
end_measure = 108

# 3. As requested, output the numbers that will form the final answer.
print(f"The opus number is: {opus_number}")
print(f"The start measure of the joke is: {start_measure}")
print(f"The end measure of the joke is: {end_measure}")

# 4. Assemble the final string in the specified format.
# The format is: composer surname, opus number, measure-range
measure_range = f"{start_measure}-{end_measure}"
final_answer = f"{composer_surname}, {opus_number}, {measure_range}"

# 5. Print the fully formatted final answer.
print("\nThe final formatted answer is:")
print(final_answer)