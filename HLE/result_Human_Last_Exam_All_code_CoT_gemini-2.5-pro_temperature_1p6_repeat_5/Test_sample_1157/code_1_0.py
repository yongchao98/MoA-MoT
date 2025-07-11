def generate_management_plan():
    """
    This function outlines the management plan for the clinical case provided.
    It prints the plan step-by-step, incorporating the patient's specific data.
    """
    # Patient's Clinical Data
    hba1c = 7.5
    time_post_accident_hours = 28
    lost_teeth = "Left Lateral Incisor, Canine, and First Premolar"
    snb_angle_degrees = 88

    print("### CLINICAL MANAGEMENT PLAN ###")

    # Part 1: Immediate and Urgent Management
    print("\n--- 1. IMMEDIATE & URGENT MANAGEMENT ---")
    print("Priority must be given to the patient's systemic health before definitive dental treatment.")
    print("\n**A. Systemic Medical Management (Highest Priority):**")
    print(f"- The patient's glycated hemoglobin (HbA1c) of {hba1c}% indicates uncontrolled diabetes mellitus. This condition severely compromises wound healing and increases the risk of post-traumatic infection.")
    print("- An immediate referral to a physician or endocrinologist is mandatory to stabilize her condition. Management will involve initiating diabetes treatment (e.g., insulin, oral hypoglycemics), regular blood glucose monitoring, and patient education.")
    
    print("\n**B. Urgent Dental Management (after medical assessment):**")
    print(f"- Given the {time_post_accident_hours}-hour delay since the accident, replantation of the avulsed teeth is not a viable option.")
    print("- **Wound Care:** The primary goal is to prevent infection. This includes gentle debridement and irrigation of the tooth sockets, assessment for alveolar bone fractures via radiographs (Panoramic or CBCT), and suturing of soft tissue lacerations.")
    print("- **Infection Control:** Systemic antibiotics are crucial due to the delay, wound contamination, and the patient's diabetic status. A tetanus booster should also be administered if required.")
    
    # Part 2: Cells of Interest
    print("\n--- 2. CELLS OF INTEREST ---")
    print("The key cells involved in her healing process are:")
    print("- **Immune Cells (Neutrophils and Macrophages):** Essential for fighting infection. Their function is impaired by high blood sugar.")
    print("- **Fibroblasts:** Responsible for producing collagen for soft tissue healing. Their function is also hindered by uncontrolled diabetes.")
    print("- **Osteoblasts and Osteoclasts:** These cells are critical for healing and remodeling the bone in the tooth sockets. Diabetes can negatively affect this process.")

    # Part 3: Prosthodontic (Denture) Replacement Plan
    print("\n--- 3. LONG-TERM DENTURE REPLACEMENT PLAN ---")
    print("This phase should only begin after her diabetes is well-controlled (target HbA1c < 7.0%).")
    print(f"The patient's SNB angle of {snb_angle_degrees}° indicates a prognathic mandible (Skeletal Class III), which must be factored into the denture design.")

    print("\n**A. Type of Denture and Material:**")
    print("- **Denture Type:** A **Removable Partial Denture (RPD)** is the most suitable initial choice. It is conservative, facilitates good oral hygiene (vital for diabetics), and is cost-effective.")
    print("- **Material:** A **Cobalt-Chromium cast metal framework** is recommended for its strength, rigidity, and thin profile, combined with an **acrylic resin** base and denture teeth.")

    print("\n**B. Abutment Teeth and Rationale:**")
    print(f"- **Lost Teeth:** {lost_teeth}.")
    print("- **Abutment Teeth:** The primary abutments would be the teeth adjacent to the space: the **Left Second Premolar** and the **Left Central Incisor**. For stability, abutments on the opposite side (e.g., Right Canine/Premolar) should be included for **cross-arch stabilization**.")
    print("- **Rationale:** This design distributes biting forces across multiple strong teeth, preventing overloading. Cross-arch stabilization is crucial to resist dislodging forces, especially considering her Skeletal Class III jaw relationship. The health of these abutment teeth must be confirmed before proceeding.")
    
    print("\n--- FINAL NUMERICAL SUMMARY ---")
    print("The management plan is based on the following key clinical numbers:")
    print(f"Equation of Factors: Management(HbA1c = {hba1c}%, Delay = {time_post_accident_hours} hours, SNB = {snb_angle_degrees}°)")


if __name__ == '__main__':
    generate_management_plan()
