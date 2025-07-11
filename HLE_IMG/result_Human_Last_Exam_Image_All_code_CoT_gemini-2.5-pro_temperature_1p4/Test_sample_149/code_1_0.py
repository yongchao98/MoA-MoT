def provide_diagnosis():
    """
    This script details the logical steps to reach a diagnosis based on the provided
    clinical and histopathological data.
    """
    # Step 1: Analyze Clinical Presentation
    patient_age = 28
    location = "sella turcica"
    print("Step 1: Clinical Analysis")
    print(f"A {patient_age}-year-old female presents with a mass in the {location}.")
    print("This location is near the sphenoid sinus, a site prone to fungal infections that can extend intracranially.")
    print("-" * 30)

    # Step 2: Analyze Histopathology
    stain = "Periodic acid-Schiff (PAS)"
    hyphae_type = "septate (with cross-walls)"
    branching_angle = 45
    print("Step 2: Histopathological Image Analysis")
    print(f"The image is stained with a {stain}-like stain, highlighting fungal elements in magenta.")
    print("Key microscopic findings include:")
    print(f" - Filamentous fungal hyphae that are {hyphae_type}.")
    print(f" - The hyphae exhibit dichotomous branching at acute angles (approximately {branching_angle} degrees).")
    print(" - A dense inflammatory infiltrate is present in the surrounding tissue.")
    print("-" * 30)

    # Step 3: Conclusion
    diagnosis = "Aspergillosis"
    print("Step 3: Conclusion")
    print("The morphology of acutely branching, septate hyphae is the hallmark of Aspergillus species.")
    print("The clinical presentation of a sellar mass combined with these histopathological features confirms the diagnosis.")
    print("\nFinal Diagnosis:")
    print(diagnosis)

if __name__ == "__main__":
    provide_diagnosis()