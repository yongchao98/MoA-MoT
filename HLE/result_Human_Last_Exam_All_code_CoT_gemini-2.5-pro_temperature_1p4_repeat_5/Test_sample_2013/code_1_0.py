# This script will identify and print the details of a famous musical joke in a piano concerto.
# The joke involves a cello flourish that becomes desynchronized and then recovers.

# 1. Define the composer's surname.
# The piece is the Piano Concerto No. 5, "Emperor", by Ludwig van Beethoven.
composer_surname = "Beethoven"

# 2. Define the opus number of the concerto.
# Piano Concerto No. 5 is Opus 73.
opus_number = 73

# 3. Define the measure range for the musical joke.
# The cello flourish begins in m. 355.
# The cello is off-beat for two measures (mm. 357-358).
# The cello skips a beat and recovers in m. 359.
# The full range spanning the event is from the start of the flourish to the recovery.
start_measure = 355
end_measure = 359

# 4. Format the final answer string with the three required elements.
# The format is: composer surname, opus number, range of measures.
final_answer = f"{composer_surname}, {opus_number}, {start_measure}-{end_measure}"

# 5. Print the final formatted answer.
print(final_answer)