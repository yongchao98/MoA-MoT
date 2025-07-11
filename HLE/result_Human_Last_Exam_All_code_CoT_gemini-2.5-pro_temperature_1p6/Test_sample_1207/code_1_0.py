import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = buffer = io.StringIO()

def solve_medical_case():
    """
    Analyzes the clinical case and identifies the most likely imaging finding.
    """
    # Patient data from the prompt
    age = 44
    gravida = 2
    para = 1

    # Rationale for the decision
    explanation = (
        "The patient presents with a multi-system inflammatory disorder affecting the nervous system (transient vision loss, headaches, hearing loss), lungs (dyspnea), joints (pain), and skin (painful leg area). "
        "This clinical picture is highly suggestive of Sarcoidosis.\n\n"
        "Neurosarcoidosis, a common manifestation, can cause inflammation of the brain's covering (leptomeninges) and cranial nerves, explaining the headaches, vision, and hearing issues. "
        "Therefore, an MRI of the brain is the most appropriate imaging modality to investigate these specific symptoms.\n\n"
        "The finding of 'Leptomeningeal enhancement with \"snowball\" hyperintensities' is a classic and specific sign of neurosarcoidosis on an MRI scan, making it the most expected finding."
    )
    
    print(explanation)
    
    # Printing the numbers from the case as requested by the prompt instructions
    print("\nPatient Demographic Numbers:")
    print(f"Age: {age}")
    print(f"Gravida: {gravida}")
    print(f"Para: {para}")
    
    # The final answer in the required format
    final_answer = "<<<B>>>"
    print(f"\nFinal Answer: {final_answer}")

# Run the solver
solve_medical_case()

# Restore original stdout and get the captured output
sys.stdout = original_stdout
output = buffer.getvalue()

# Print the captured output to the actual console
print(output)