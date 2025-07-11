def explain_bulk_experiment_limitation():
    """
    This function explains the fundamental limitation of bulk calorimetric experiments
    in the context of nucleic acid thermodynamics.
    """
    title = "Analysis of Limitations in Bulk Calorimetry"
    print(title)
    print("=" * len(title))

    explanation = """
The question concerns a fundamental limitation of studying nucleic acid (NA) thermodynamics using bulk melting experiments. A "bulk" experiment measures the average signal from a massive population of molecules.

Let's evaluate the given choices:

*   A. Heat capacity change is assumed to be zero.
*   B. The NNBP parameters are T-independent.
    - These are features of the simplest *analytical models* (like the van't Hoff approximation) used to interpret the data, not inherent limitations of the *experimental method* itself. More sophisticated models can account for these factors if the data is precise enough.

*   D. Temperature oscillations in bulk calorimetry are too large to capture T-dependence.
*   E. Temperature cannot be controlled in calorimetric experiments.
    - These statements are incorrect. Modern calorimeters offer very precise temperature control and measurement.

*   C. Impossibility to capture heterogeneity in bulk experiments.
    - This is the core, fundamental limitation. A bulk measurement provides an ensemble average. It cannot distinguish between different molecular states or conformations that may co-exist in the sample (e.g., correctly folded helices, misfolded structures, alternative secondary structures, or aggregates). All these different sub-populations contribute to the final, averaged melting curve, masking the underlying complexity and heterogeneity. This is a limitation that persists even with perfect experimental technique.
"""
    print(explanation)

explain_bulk_experiment_limitation()
<<<C>>>