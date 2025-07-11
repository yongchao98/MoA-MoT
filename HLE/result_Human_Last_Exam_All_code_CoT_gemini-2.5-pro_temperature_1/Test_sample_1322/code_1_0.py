def predict_gene_persistence(selection_coefficient, effective_population_size):
    """
    Models the persistence of a gene based on the balance between natural selection and genetic drift.

    In the context of genomic decay, a gene's persistence depends on whether the
    selective pressure to keep it is strong enough to overcome random genetic drift.
    A common approximation for this threshold is when the selection coefficient (s)
    is greater than the reciprocal of the effective population size (1/Ne). This comparison
    is a proxy for the efficiency of natural selection.

    Args:
        selection_coefficient (float): The fitness advantage (s) of retaining the gene.
        effective_population_size (int): The effective population size (Ne).
    """
    print(f"--- Analyzing Gene Persistence ---")
    print(f"Scenario with Selection Coefficient (s) = {selection_coefficient} and Effective Population Size (Ne) = {effective_population_size}")

    # The efficiency of selection is challenged by genetic drift.
    # The strength of drift is inversely proportional to population size.
    strength_of_drift = 1 / effective_population_size

    # The "equation" we are evaluating is whether selection is stronger than drift.
    # We will print each component of this comparison.
    print(f"The strength of selection is represented by s = {selection_coefficient}")
    print(f"The strength of genetic drift is approximated by 1/Ne = 1/{effective_population_size} = {strength_of_drift:.6f}")

    # If s > 1/Ne, selection is considered "efficient" enough to act on this gene.
    # If s <= 1/Ne, selection is "inefficient" and drift dominates.
    if selection_coefficient > strength_of_drift:
        print(f"\nResult: Since s ({selection_coefficient}) > 1/Ne ({strength_of_drift:.6f}), selection is relatively efficient.")
        print("Conclusion: The genomic fragment is likely to PERSIST.")
    else:
        print(f"\nResult: Since s ({selection_coefficient}) <= 1/Ne ({strength_of_drift:.6f}), genetic drift is likely to overpower selection.")
        print("Conclusion: The genomic fragment is likely to be LOST during genomic decay.")

# --- Example Scenarios ---

# Case 1: A weakly beneficial gene in a large population.
# Selection is efficient enough to preserve the gene.
print("Case 1: Large Population (High Selection Efficiency)")
predict_gene_persistence(selection_coefficient=0.001, effective_population_size=50000)

print("\n" + "="*40 + "\n")

# Case 2: The same weakly beneficial gene in a small endosymbiotic population.
# Drift is much stronger, making selection inefficient for this gene.
print("Case 2: Small Population (Low Selection Efficiency)")
predict_gene_persistence(selection_coefficient=0.001, effective_population_size=500)