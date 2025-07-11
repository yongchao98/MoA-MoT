def triage_spine_patients():
    """
    This function prioritizes patients based on their spinal injuries
    and prints the reasoning and the final order.
    """
    
    # Patient data
    patients = {
        "Patient 1": "Severe burst fracture of L2, no neurologic deficits.",
        "Patient 2": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.",
        "Patient 3": "Split fracture of L2, with mildly disordered pelvic functions."
    }

    print("Spine Surgery Triage Rationale:")
    print("---------------------------------")
    
    # Priority 1: Neurologic Deficit
    print("1. Highest Priority: Patient 3")
    print(f"   Diagnosis: {patients['Patient 3']}")
    print("   Reasoning: The presence of 'disordered pelvic functions' indicates an active neurologic deficit (cauda equina syndrome). This is a surgical emergency requiring immediate decompression to prevent permanent nerve damage.")
    
    # Priority 2: High Instability
    print("\n2. Second Priority: Patient 1")
    print(f"   Diagnosis: {patients['Patient 1']}")
    print("   Reasoning: A 'severe burst fracture' is a highly unstable injury. Even without current neurologic deficits, there is a significant risk of future spinal cord or nerve root injury due to mechanical instability. Surgery is urgent to stabilize the spine.")

    # Priority 3: Lower Instability
    print("\n3. Lowest Priority: Patient 2")
    print(f"   Diagnosis: {patients['Patient 2']}")
    print("   Reasoning: A compression fracture with 'mild' spondylolisthesis is the least unstable of the three conditions. With no neurologic deficits, the risk of imminent deterioration is lower, making this case less urgent than the others.")

    # Final Prioritization Order
    final_order = ["Patient 3", "Patient 1", "Patient 2"]
    
    print("\n--- Final Prioritization Order ---")
    print(f"The final surgical priority from highest to lowest is:")
    print(f"{final_order[0]} > {final_order[1]} > {final_order[2]}")

# Execute the function
triage_spine_patients()
<<<F>>>