import textwrap

def solve_clinical_case():
    """
    Analyzes a clinical case to identify the most likely chromosomal abnormality.
    """
    patient_features = {
        "Craniofacial": ["Cleft palate", "Microcephaly", "Frontal bossing", "Prominent eyes", "Midface hypoplasia", "Micrognathia", "Dysplastic ears"],
        "Skeletal": ["Broadened shoulders", "Short stature", "Clinodactyly"],
        "Growth/Development": ["Low BMI for age", "Delayed developmental milestones", "Intellectual disability"],
        "Psychiatric/Speech": ["Incomprehensible speech"],
        "History": ["Born preterm"]
    }

    # Analysis of options
    print("Step 1: Analyzing the patient's key features.")
    key_symptoms = ["Cleft palate", "Developmental delay/Intellectual Disability", "Midface hypoplasia", "Micrognathia", "Dysplastic ears", "Speech impairment"]
    print(f"The most distinguishing features are: {', '.join(key_symptoms)}.\n")

    print("Step 2: Evaluating the chromosomal abnormality options.")

    # Option A: Chromosome 3
    print("A. Chromosome 3: Associated with conditions like 3q29 microdeletion syndrome. While it causes developmental delay and microcephaly, the specific combination of cleft palate and facial features in the patient is less characteristic.")

    # Option B: Chromosome 22
    analysis_b = """
B. Chromosome 22: Strongly associated with 22q11.2 deletion syndrome (DiGeorge/Velocardiofacial syndrome).
   - Cleft palate: A classic feature. (Match)
   - Intellectual disability and developmental delay: Common. (Match)
   - Facial features (midface hypoplasia, micrognathia, dysplastic ears): Hallmarks of the syndrome. (Match)
   - Speech issues: Very common due to palatal defects and developmental delay. (Match)
   - Short stature and clinodactyly can also be associated. (Match)
This is a very strong match for the patient's presentation.
"""
    print(textwrap.dedent(analysis_b))

    # Option C: Chromosome 21
    print("C. Chromosome 21: Associated with Trisomy 21 (Down Syndrome). While there is an intellectual disability, the specific facial dysmorphisms (flattened profile, up-slanting eyes) are different from the patient's. This is not a good match.")

    # Option D: Chromosome 2
    print("D. Chromosome 2: Associated with various abnormalities like 2q37 deletion syndrome. While there can be intellectual disability, the constellation of symptoms, especially the prominent palatal and facial features, is not the most classic presentation.")

    # Option E: Chromosome 13
    print("E. Chromosome 13: Associated with Trisomy 13 (Patau Syndrome). This is a severe condition, and affected individuals rarely survive past infancy. The patient is 15 years old, making this diagnosis highly unlikely.\n")

    print("Step 3: Conclusion.")
    print("Based on the analysis, the patient's combination of cleft palate, characteristic facial features, and developmental/speech delays is most consistent with a disorder linked to Chromosome 22, specifically 22q11.2 deletion syndrome.")
    print("\nTherefore, the expected chromosomal abnormality is on Chromosome 22.")
    print("\nFinal Answer: B")


solve_clinical_case()