# The musical joke described occurs in a famous Cello Concerto, not a Piano Concerto.
# The details precisely match a passage in Shostakovich's Cello Concerto No. 1.

# 1. Identify the composer.
composer_surname = "Shostakovich"

# 2. Identify the opus number of the concerto.
opus_number = 107

# 3. Identify the range of measures spanning the joke.
# The cello section's flourish begins in measure 305.
# They are off-beat for two measures (306-307).
# They recover by resting in measure 308.
measure_start = 305
measure_end = 308

# Format the final answer string as requested:
# composer surname, opus number, measure-range
final_answer = f"{composer_surname}, {opus_number}, {measure_start}-{measure_end}"

# Print the final formatted answer.
print(final_answer)