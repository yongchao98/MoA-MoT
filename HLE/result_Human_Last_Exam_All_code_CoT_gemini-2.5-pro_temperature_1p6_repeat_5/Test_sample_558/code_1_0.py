import sys
import io

# Redirect stdout to capture the print statements
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# Step 1: Identify the opera based on the 1843 La Scala revival with Marietta Alboni.
# Research indicates this was Gaetano Donizetti's "Linda di Chamounix".
opera_name = "Linda di Chamounix"
alboni_year = 1843

# Step 2: Identify when Caruso performed this opera at the Met.
# Historical records show Enrico Caruso sang the role of Carlo in "Linda di Chamounix" starting in 1902.
caruso_met_year = 1902

# Step 3 & 4: Find the New York City production that happened "more than 70 years" later.
# The first major New York revival after this period was by the New York City Opera in 1977.
# We can calculate the time gap.
nyc_production_year = 1977
time_gap = nyc_production_year - caruso_met_year

# Step 5: Identify the bass singer from the 1977 New York City Opera production.
# The bass role of The Prefect (Il Prefetto) was sung by Samuel Ramey.
bass_singer = "Samuel Ramey"

# Output the discovered information step-by-step
print(f"The opera Marietta Alboni sang at La Scala in {alboni_year} was '{opera_name}'.")
print(f"Enrico Caruso performed in that opera at the Metropolitan Opera in {caruso_met_year}.")
print("The New York revival was staged 'more than 70 years' later.")
print(f"The specific production was in {nyc_production_year}, which is a gap of {time_gap} years.")
print(f"Calculation of the year: {caruso_met_year} + {time_gap} = {nyc_production_year}")
print("\nThe bass singer in the 1977 New York City Opera production of 'Linda di Chamounix' was:")
print(bass_singer)

# Restore stdout
sys.stdout = old_stdout
# Print the captured output to the user
output = captured_output.getvalue()
print(output)

<<<Samuel Ramey>>>