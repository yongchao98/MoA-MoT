# 1. Define the components of the answer based on musical analysis.
# Composer: The piece is by Ludwig van Beethoven.
composer_surname = "Beethoven"

# Opus Number: The piano concerto is his No. 2 in B-flat major, which is Opus 19.
opus_number = 19

# Measure Range: In the 3rd movement (Rondo), the humorous cello passage begins its flourish
# in measure 286, is off-beat in measures 287-288, and recovers in measure 289.
# The requested range spans from the beginning of the flourish to the recovery (inclusive).
start_measure = 286
end_measure = 289

# 2. Format the final string as a comma-separated list.
# The measure range is represented as a hyphenated pair.
measure_range_str = f"{start_measure}-{end_measure}"
final_answer = f"{composer_surname}, {opus_number}, {measure_range_str}"

# 3. Print the result.
print(final_answer)