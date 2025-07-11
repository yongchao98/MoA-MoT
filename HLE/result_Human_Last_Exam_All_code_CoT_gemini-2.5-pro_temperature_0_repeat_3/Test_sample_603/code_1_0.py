def solve_biophysics_question():
    """
    This script analyzes a multiple-choice question about the limitations of
    bulk experiments in nucleic acid thermodynamics and prints the correct answer with reasoning.
    """

    question = "What is a fundamental limitation of studying nucleic acid (NA) thermodynamics using bulk melting experiments, even under ideal experimental conditions?"

    options = {
        'A': "Heat capacity change is assumed to be zero.",
        'B': "The NNPB parameters are T-independent.",
        'C': "Impossibility to capture heterogeneity in bulk experiments.",
        'D': "Temperature oscillations in bulk calorimetry are too large to capture T-dependence.",
        'E': "Temperature cannot be controlled in calorimetric experiments."
    }

    print("Analyzing the question: " + question)
    print("-" * 30)
    print("Step-by-step analysis of the options:")
    print("\n1. Options A and B relate to the thermodynamic model, not the experiment itself.")
    print("   - The assumption of zero heat capacity change (ΔCp = 0) and temperature-independent parameters (ΔH, ΔS) is a feature of the simplest two-state model.")
    print("   - It's a modeling choice, not an inherent limitation of the bulk measurement technique. More complex models can be applied to the same bulk data.")

    print("\n2. Options D and E are factually incorrect regarding modern instrumentation.")
    print("   - Modern calorimeters and spectrophotometers provide extremely precise temperature control and ramping. Without this control, the experiment would be impossible.")
    print("   - Therefore, these are not fundamental limitations.")

    print("\n3. Option C describes a core principle of bulk vs. single-molecule measurements.")
    print("   - A 'bulk' experiment measures the ensemble average of a massive population of molecules.")
    print("   - This average masks any heterogeneity within the sample. For example, if a small fraction of molecules is misfolded or populates an intermediate state, the bulk measurement will only reflect the average behavior, not the distinct subpopulations.")
    print("   - This inability to resolve molecular-level differences is a fundamental limitation of any bulk technique.")

    print("-" * 30)
    print("Conclusion: The most fundamental limitation inherent to the 'bulk' nature of the experiment is its inability to capture molecular heterogeneity.")
    print(f"Therefore, the correct answer is C: {options['C']}")

solve_biophysics_question()