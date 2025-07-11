def diagnose_from_findings(findings):
    """
    This function simulates the diagnostic logic based on histological features.
    """
    diagnosis = "Undetermined"

    # Key features for diagnosing invasive fungal infections
    has_fungal_hyphae = findings.get("fungal_hyphae_present", False)
    hyphae_width = findings.get("hyphae_width")
    hyphae_septa = findings.get("hyphae_septa")
    branching_angle = findings.get("branching_angle")
    location = findings.get("location")

    print("Analyzing Clinical and Histological Findings:")
    print(f"- Location of mass: {location}")
    print(f"- Fungal hyphae present: {has_fungal_hyphae}")
    print(f"- Hyphae width: {hyphae_width}")
    print(f"- Hyphae septation: {hyphae_septa}")
    print(f"- Hyphae branching angle: Approximately {branching_angle} degrees")
    print("-" * 20)

    # Logic to differentiate common invasive fungal infections
    if has_fungal_hyphae:
        if hyphae_width == "broad" and hyphae_septa == "pauciseptate" and branching_angle == 90:
            diagnosis = "Mucormycosis"
        elif hyphae_width == "narrow" and hyphae_septa == "septate" and branching_angle == 45:
            diagnosis = "Aspergillosis"
        else:
            diagnosis = "Unspecified Fungal Infection"
    else:
        # If no fungi were seen, other diagnoses for a sellar mass would be considered.
        if location == "sella turcica":
            diagnosis = "Differential includes Pituitary Adenoma, but this contradicts the image."
        else:
            diagnosis = "No fungal elements identified."
            
    return diagnosis

# Patient data extracted from the prompt and image analysis
patient_findings = {
    "location": "sella turcica",
    "fungal_hyphae_present": True,
    "hyphae_width": "broad",
    "hyphae_septa": "pauciseptate",  # also known as aseptate or sparsely septate
    "branching_angle": 90
}

# Run the diagnosis
final_diagnosis = diagnose_from_findings(patient_findings)

print(f"Final Diagnosis: {final_diagnosis}")
