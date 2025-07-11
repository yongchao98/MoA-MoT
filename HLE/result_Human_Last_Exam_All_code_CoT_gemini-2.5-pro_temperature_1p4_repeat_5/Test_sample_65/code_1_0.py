def triage_spine_patients():
    """
    Prioritizes three spine injury patients based on clinical descriptions.
    The function explains the reasoning and prints the final priority order.
    """
    
    # Patient data
    patients = {
        'Patient 1': 'Severe burst fracture of L2, no neurologic deficits.',
        'Patient 2': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.',
        'Patient 3': 'Split fracture of L2, with mildly disordered pelvic functions.'
    }
    
    # Priority is determined by:
    # 1. Presence of neurologic deficit (highest priority)
    # 2. Degree of spinal instability
    
    print("Evaluating patients based on surgical indications:\n")
    
    # Highest Priority: Patient 3
    print("Priority 1: Patient 3")
    print("Reasoning: This patient shows signs of active neurologic deficit ('mildly disordered pelvic functions'), suggesting cauda equina syndrome. This is a surgical emergency to prevent permanent neurological damage.")
    print("-" * 50)

    # Second Priority: Patient 1
    print("Priority 2: Patient 1")
    print("Reasoning: This patient has a 'severe burst fracture,' which signifies major spinal instability. While there are no current deficits, the risk of future neurologic injury or deformity is very high, making surgical stabilization urgent.")
    print("-" * 50)
    
    # Third Priority: Patient 2
    print("Priority 3: Patient 2")
    print("Reasoning: This patient has a compression fracture with 'mild' spondylolisthesis. This represents a less severe and less unstable injury pattern compared to a severe burst fracture. With no neurologic deficits, this case is the least urgent of the three.")
    print("-" * 50)

    # Final prioritized list
    final_order = ["Patient 3", "Patient 1", "Patient 2"]
    
    print("Final Triage Order (from top priority to lowest priority):")
    print(f"{final_order[0]}, {final_order[1]}, {final_order[2]}")

# Run the triage
triage_spine_patients()
<<<F>>>