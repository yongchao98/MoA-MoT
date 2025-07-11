# The musical joke in question is from Beethoven's Piano Concerto No. 5, Op. 73.
# This event occurs in the third movement (Rondo).

# 1. Define the composer's surname.
composer_surname = "Beethoven"

# 2. Define the opus number of the concerto.
opus_number = 73

# 3. Define the start and end measures of the joke.
#    - The flourish begins in measure 343.
#    - The recovery, by skipping a beat, happens in measure 346.
start_measure = 343
end_measure = 346

# 4. Format the measure range as a hyphenated string.
measure_range = f"{start_measure}-{end_measure}"

# 5. Combine all elements into the final comma-separated string and print it.
#    The output will be the composer, opus number, and the full measure range of the joke.
print(f"{composer_surname}, {opus_number}, {measure_range}")