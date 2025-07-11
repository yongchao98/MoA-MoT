# This script formats and prints the details of a famous musical joke in a classical piano concerto.
# The information identifies the composer, the work, and the specific measures where the event occurs.

composer_surname = "Haydn"
concerto_number = 11 # From the Hoboken-Verzeichnis catalog number: Hob. XVIII:11
joke_start_measure = 136 # Cello/bassoon play a 7-note phrase (1 beat too long)
recovery_measure = 139 # Cello/bassoon play a 5-note phrase to get back on beat

# The range spans from the start of the flourish to the measure where the beat is recovered.
measure_range = f"{joke_start_measure}-{recovery_measure}"

# Print the final formatted answer.
print(f"{composer_surname}, {concerto_number}, {measure_range}")