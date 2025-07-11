# The musical joke described is a famous passage in a concerto by composer Emil von Sauer.

# 1. Identify the composer's surname.
composer_surname = "Sauer"

# 2. Identify the work. It is his Piano Concerto No. 1. This work does not have an opus number,
# so we will use the concerto number as its identifier.
concerto_identifier = 1

# 3. Identify the measure range of the joke in the third movement.
# The cello flourish begins in measure 139.
start_measure = 139
# The cellos are off-beat for the next two measures and recover by skipping a beat in measure 143.
end_measure = 143

# Print the final formatted answer.
print(f"{composer_surname}, {concerto_identifier}, {start_measure}-{end_measure}")