import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
string_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = string_buffer

# --- Start of the student's code ---

# Specimen counts from the collection
dixidae_specimens = 100
simuliidae_specimens = 1101
belpharoceridae_specimens = 100
rhagionidae_specimens = 101
tabanidae_specimens = 201

# Number of prolegs/parapodia per larva for each family
# Dixidae larvae have 2 pairs of prolegs.
dixidae_prolegs = 4
# Simuliidae larvae have 1 anterior and 1 posterior proleg.
simuliidae_prolegs = 2
# Belpharoceridae larvae have 6 ventral suckers/prolegs.
belpharoceridae_prolegs = 6
# Rhagionidae (Vermileoninae) larvae are apodous (no prolegs).
rhagionidae_prolegs = 0
# Tabanidae larvae have prominent fleshy prolegs (creeping welts), simplified to 4 pairs.
tabanidae_prolegs = 8

# Calculate the subtotal of prolegs for each family
total_dixidae = dixidae_specimens * dixidae_prolegs
total_simuliidae = simuliidae_specimens * simuliidae_prolegs
total_belpharoceridae = belpharoceridae_specimens * belpharoceridae_prolegs
total_rhagionidae = rhagionidae_specimens * rhagionidae_prolegs
total_tabanidae = tabanidae_specimens * tabanidae_prolegs

# Calculate the grand total
grand_total = total_dixidae + total_simuliidae + total_belpharoceridae + total_rhagionidae + total_tabanidae

# Print the full equation showing the calculation for each family
print(f"({dixidae_specimens} * {dixidae_prolegs}) + ({simuliidae_specimens} * {simuliidae_prolegs}) + ({belpharoceridae_specimens} * {belpharoceridae_prolegs}) + ({rhagionidae_specimens} * {rhagionidae_prolegs}) + ({tabanidae_specimens} * {tabanidae_prolegs})")
print(f"= {total_dixidae} + {total_simuliidae} + {total_belpharoceridae} + {total_rhagionidae} + {total_tabanidae}")
print(f"= {grand_total}")

# --- End of the student's code ---

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
captured_output = string_buffer.getvalue()

# Print the captured output to the user
print(captured_output)

# Extract the final numerical answer for the grading format
final_answer = grand_total
print(f"<<<{final_answer}>>>")