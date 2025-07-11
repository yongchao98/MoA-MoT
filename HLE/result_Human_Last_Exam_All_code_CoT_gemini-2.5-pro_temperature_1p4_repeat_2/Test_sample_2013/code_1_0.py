import sys

# Step 1: Identify the composer and the work.
# The composer is Joseph Haydn.
# The piece is the Piano Concerto No. 11 in D Major, Hob. XVIII:11.
# The musical joke occurs in the 3rd movement, "Rondo all'Ungarese".
composer_surname = "Haydn"
work_id = "Hob. XVIII:11"

# Step 2: Identify the specific measures.
# In the finale, near the end, the cellos and basses play a short accompanying figure.
# In measure 147, they play the figure as a solo lead-in.
# The orchestra expects to enter on the downbeat of measure 148, but the cellos mistakenly
# play their figure again, creating a clash and making them off-beat.
# They remain off-beat through measure 149.
# In measure 150, they skip a beat (play a rest) to realign with the orchestra.
# Therefore, the entire event spans from the start of the flourish to the recovery.
start_measure = 147
end_measure = 150

# Step 3: Format the answer as a comma-separated list of three elements.
# The elements are: Composer Surname, Opus Number (or equivalent catalog number), and the measure range.
measure_range = f"{start_measure}-{end_measure}"
final_answer = f"{composer_surname}, {work_id}, {measure_range}"

# Step 4: Print the final formatted string.
print(final_answer)

# Ensure the final answer is also available for the <<<...>>> format if needed,
# by adding it to a known location for the wrapper.
# This part is just for demonstrating the variables are correct.
# In a real script, only the print statement is necessary for output.
sys.stdout.flush() # The platform will capture the printed output.