import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# The house is at 433 E 2100th Ave, which has a mailing address in Beecher, IL.
# However, Beecher's population is ~4,713 (2020 census), which is not over 5,000.
# The question includes the constraint "(The town has a population of over 5,000)".
# The nearest town in Illinois that meets this population criterion is Crete.
# Population of Crete, IL (2020 census) is 8,465.
# Therefore, assuming the question implies the nearest qualifying town, the answer is Crete.

town = "Crete"
print(town)

# Restore stdout
sys.stdout = old_stdout
output = captured_output.getvalue()

# Final print to the user
print(output)