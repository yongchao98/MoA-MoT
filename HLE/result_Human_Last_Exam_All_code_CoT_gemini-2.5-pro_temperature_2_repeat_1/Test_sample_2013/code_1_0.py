# This script identifies and formats information about a famous musical joke in a classical piano concerto.

# 1. The composer is Ludwig van Beethoven.
composer_surname = "Beethoven"

# 2. The piece is Piano Concerto No. 5, which is Opus 73.
opus_number = 73

# 3. The joke, where the cellos/basses enter a beat early and then recover,
#    begins with their entry in measure 107 and resolves in measure 111.
measure_range = "107-111"

# 4. Format the final answer as a comma-separated list as requested.
# The elements are: composer surname, opus number, and the hyphenated measure range.
final_answer = f"{composer_surname}, {opus_number}, {measure_range}"

# Print the final result to the console.
print(final_answer)