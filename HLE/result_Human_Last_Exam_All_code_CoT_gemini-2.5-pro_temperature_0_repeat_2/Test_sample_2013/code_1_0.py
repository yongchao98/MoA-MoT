# This script identifies and formats the details of a famous musical joke
# by Ludwig van Beethoven.

# 1. Identify the composer.
composer_surname = "Beethoven"

# 2. Identify the opus number of the piano concerto.
# The work is the Piano Concerto No. 1 in C major.
opus_number = 15

# 3. Identify the measure range of the joke in the third movement (Rondo).
# The cello and bass flourish begins in measure 403, creating an off-beat
# effect through measure 405. They skip a beat to recover in measure 406.
measure_range = "403-406"

# 4. Format and print the final answer as a comma-separated list.
print(f"{composer_surname}, {opus_number}, {measure_range}")