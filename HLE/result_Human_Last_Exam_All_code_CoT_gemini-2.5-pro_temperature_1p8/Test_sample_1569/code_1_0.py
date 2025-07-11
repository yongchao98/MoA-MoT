import textwrap

def solve_clinical_vignette():
    """
    This function analyzes the provided clinical case to determine the correct diagnosis.
    """

    # Step 1: Define the patient's data from the prompt.
    patient_data = {
        "Age": "27-year-old",
        "Symptoms": ["4 days of fever", "headaches", "myalgia", "disorientation", "heart murmur"],
        "History": "recent camping trip to Oklahoma",
        "Labs": "elevated IgM with negative IgG Lyme serology titer"
    }

    # Step 2: Explain the reasoning based on the provided data.
    print("Analyzing the clinical case step-by-step:")
    print("-" * 40)
    print("1. Patient's Presentation:")
    print(f"   - The patient has systemic symptoms (fever, headache, myalgia) plus specific organ involvement: neurologic (disorientation) and cardiac (heart murmur).")
    print(f"   - The history of camping in Oklahoma points towards a tick-borne illness.")
    print("-" * 40)
    print("2. Interpreting the Lab Results:")
    print(f"   - The key finding is: '{patient_data['Labs']}'.")
    print(f"   - IgM antibodies are produced during an acute (early) infection.")
    print(f"   - IgG antibodies are produced later in the course of an infection.")
    print(f"   - An elevated IgM with a negative IgG is the classic serological profile for an acute infection.")
    print("-" * 40)
    print("3. Connecting Labs to the Disease:")
    print(f"   - The prompt explicitly states it is a 'Lyme serology titer'.")
    print(f"   - Lyme disease is caused by the bacterium Borrelia burgdorferi.")
    print(f"   - Therefore, the patient has a positive IgM titer for Borrelia burgdorferi.")
    print("-" * 40)
    print("4. Conclusion:")
    print("   - The patient's clinical picture of neurologic and cardiac symptoms is characteristic of early disseminated Lyme disease.")
    print("   - The lab results confirm an acute infection with the agent that causes Lyme disease.")
    print("   - The positive titer is for Borrelia burgdorferi.")
    print("-" * 40)

    # There is no equation, so we will print the final answer choice.
    final_answer = "C. Borrelia burgdorferi"
    print(f"The correct choice is: {final_answer}")


# Execute the analysis
solve_clinical_vignette()