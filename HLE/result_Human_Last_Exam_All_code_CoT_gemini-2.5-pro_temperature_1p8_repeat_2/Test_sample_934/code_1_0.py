import sys

def analyze_esophageal_case():
    """
    Analyzes a clinical case to determine the most likely diagnosis
    by scoring risk factors and clinical findings.
    """
    # Step 1: Define Patient Data from the case vignette
    # Risk Factors
    heavy_smoking = 20  # Strong risk factor for SCC
    heavy_alcohol = 20  # Strong risk factor for SCC
    # Clinical Findings
    symptoms = 10  # Chest pain/odynophagia are common to many
    imaging_wall_thickening = 30  # Highly suggestive of an infiltrative process like cancer
    endoscopy_normal_mucosa = 0  # This is the key finding. It rules out many surface-level pathologies.

    print("Analyzing patient data based on a scoring system...\n")

    # Step 2 & 3: Score each diagnosis
    
    # A. Streptococcal esophagitis
    # Reasoning: Rare, usually has ulcers/plaques on endoscopy.
    score_strep = symptoms - 40 # Strong negative for normal endoscopy
    print("Diagnosis: Streptococcal esophagitis")
    print(f"Calculation: Symptoms ({symptoms}) - Incompatible Endoscopy (40) = {score_strep}")
    print("Reasoning: Infectious esophagitis is extremely unlikely with a normal-appearing mucosa on endoscopy.\n")

    # B. Esophageal adenocarcinoma
    # Reasoning: Main risk is GERD/Barrett's, not mentioned. While smoking is a risk, it's weaker than for SCC.
    score_adeno = (heavy_smoking / 2) + symptoms - 20 # Moderate risk from smoking, but unlikely without endoscopic evidence of Barrett's or a tumor mass.
    print("Diagnosis: Esophageal adenocarcinoma")
    print(f"Calculation: Smoking Risk ({heavy_smoking}/2) + Symptoms ({symptoms}) - No Endoscopic Mass (20) = {score_adeno}")
    print("Reasoning: While possible, the primary risk factors (chronic GERD) are absent, and a normal endoscopy makes a visible tumor less likely.\n")

    # C. Esophageal squamous cell carcinoma (SCC)
    # Reasoning: Patient has the two strongest risk factors (smoking, alcohol). Imaging fits perfectly.
    # Normal endoscopy is a classic presentation for an *infiltrative* or *submucosal* SCC.
    score_scc = heavy_smoking + heavy_alcohol + symptoms + imaging_wall_thickening
    print("Diagnosis: Esophageal squamous cell carcinoma (SCC)")
    print(f"Calculation: Heavy Smoking ({heavy_smoking}) + Heavy Alcohol ({heavy_alcohol}) + Symptoms ({symptoms}) + Imaging Finding ({imaging_wall_thickening}) = {score_scc}")
    print("Reasoning: The combination of major risk factors (smoking, alcohol) and imaging findings (wall thickening) is a classic presentation. A normal endoscopy is known to occur with infiltrative SCC that grows within the esophageal wall.\n")

    # D. GERD
    # Reasoning: GERD does not explain the significant wall thickening on imaging.
    score_gerd = symptoms - 30 # Significant negative for imaging findings atypical for GERD
    print("Diagnosis: GERD")
    print(f"Calculation: Symptoms ({symptoms}) - Incompatible Imaging (30) = {score_gerd}")
    print("Reasoning: Simple GERD would not cause significant esophageal wall thickening and luminal narrowing.\n")

    # E. Herpes esophagitis
    # Reasoning: Viral esophagitis classically shows 'volcano-like' ulcers on endoscopy.
    score_herpes = symptoms - 40 # Strong negative for normal endoscopy
    print("Diagnosis: Herpes esophagitis")
    print(f"Calculation: Symptoms ({symptoms}) - Incompatible Endoscopy (40) = {score_herpes}")
    print("Reasoning: Like other infectious causes, Herpes esophagitis is ruled out by the normal endoscopy.\n")
    
    # Step 4: Determine the highest score
    diagnoses = {
        'Streptococcal esophagitis': score_strep,
        'Esophageal adenocarcinoma': score_adeno,
        'Esophageal squamous cell carcinoma': score_scc,
        'GERD': score_gerd,
        'Herpes esophagitis': score_herpes
    }
    
    most_likely_diagnosis = max(diagnoses, key=diagnoses.get)
    
    # Step 5: Final Conclusion
    print("---CONCLUSION---")
    print(f"The most likely diagnosis is '{most_likely_diagnosis}' with a score of {diagnoses[most_likely_diagnosis]}.")
    print("This diagnosis best explains the patient's strong risk factors (heavy smoking and alcohol) and the clinical picture of an infiltrative process (wall thickening on imaging with normal surface mucosa on endoscopy).")

if __name__ == '__main__':
    analyze_esophageal_case()