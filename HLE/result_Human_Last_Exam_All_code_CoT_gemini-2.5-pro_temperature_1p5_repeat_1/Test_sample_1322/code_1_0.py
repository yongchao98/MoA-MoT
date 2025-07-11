# The Drift-Barrier Hypothesis: Modeling the persistence of genomic fragments.

# Define the variables for our model.
# Ne: Effective population size. Endosymbionts and parasites often have small population sizes.
Ne = 5000
# s: The selection coefficient. This represents the small metabolic cost (a disadvantage)
# of replicating a useless DNA fragment. We represent it as a positive number for the magnitude.
s = 0.00004

# According to the drift-barrier hypothesis, natural selection is inefficient
# when the magnitude of the selection coefficient 's' is less than the
# strength of genetic drift, which is approximated by '1/Ne'.
drift_barrier = 1 / Ne

# Print the context and the values.
print("--- Evaluating Fragment Persistence using the Drift-Barrier Hypothesis ---")
print(f"The selective disadvantage of the fragment (s): {s}")
print(f"The effective population size (Ne): {Ne}")
# This is the 'equation' part, showing the numbers that are compared.
print(f"The strength of genetic drift (1/Ne): 1 / {Ne} = {drift_barrier}")
print("-----------------------------------------------------------------------")

# Compare 's' to the drift barrier to determine the fate of the fragment.
if s < drift_barrier:
    # This is the condition where the fragment persists.
    print(f"\nResult: Since s ({s}) < 1/Ne ({drift_barrier}), natural selection is inefficient.")
    print("The weak disadvantage of the fragment is swamped by the effects of random genetic drift.")
    print("Conclusion: The genomic fragment is likely to persist in the genome.")
else:
    # This would be the case in a very large population.
    print(f"\nResult: Since s ({s}) >= 1/Ne ({drift_barrier}), natural selection is efficient.")
    print("The disadvantage of the fragment is strong enough for selection to act upon it.")
    print("Conclusion: The genomic fragment is likely to be deleted from the genome over time.")
