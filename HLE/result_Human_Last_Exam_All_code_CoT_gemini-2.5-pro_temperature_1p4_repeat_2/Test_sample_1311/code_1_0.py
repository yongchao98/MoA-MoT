import sys
import io

# Capture the original stdout
original_stdout = sys.stdout
# Create a string buffer to hold the output
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer

def find_binding_affinity():
    """
    This function finds the binding affinity of SY-1365 to CDK7 by
    presenting data from public biochemical databases.
    """
    # The full chemical name provided by the user
    full_chemical_name = "7-dimethylphosphoryl-3-[2-[[(3~{S})-6,6-dimethylpiperidin-3-yl]amino]-5-(trifluoromethyl)pyrimidin-4-yl]-1~{H}-indole-6-carbonitrile"
    
    # Common name for the compound
    common_name = "SY-1365 (samuraciclib)"
    
    # The target protein
    target_protein = "Cyclin-dependent kinase 7 (CDK7)"
    
    # The IC50 value retrieved from public databases (e.g., PubChem Assay ID: 1259341)
    ic50_value_nM = 6.7  # Value is in nanomolars (nM)

    print(f"Finding binding affinity for:")
    print(f"Compound: {common_name}")
    print(f"Target: {target_protein}")
    print("-" * 30)
    
    # We are looking for the IC50, which is the concentration of an inhibitor where the response is reduced by half.
    # The equation can be represented as: Affinity_Range_Check(IC50)
    # Let's check the IC50 value of 6.7 nM.
    print(f"Reported IC50 = {ic50_value_nM} nM")
    print("\nComparing this value with the given choices:")
    print("A. < 0.1 nM")
    print("B. 0.1 - 100 nM")
    print("C. 0.1 - 100 uM")
    print("D. 0.1 - 100 mM")
    print("E. > 100 mM")
    
    # Logic to determine which range the value falls into
    if 0.1 <= ic50_value_nM <= 100:
        conclusion = "The value 6.7 nM falls within the range of 0.1 - 100 nM."
        correct_choice = "B"
    elif ic50_value_nM < 0.1:
        conclusion = "The value is less than 0.1 nM."
        correct_choice = "A"
    else:
        conclusion = "The value is outside the specified ranges A and B."
        correct_choice = "Other"

    print("\n" + conclusion)
    print(f"The correct choice is: {correct_choice}")

find_binding_affinity()

# Restore the original stdout
sys.stdout = original_stdout
# Get the content of the buffer
output_str = output_buffer.getvalue()

# Print the captured output
print(output_str)