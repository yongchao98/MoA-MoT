import sys
import io

# Helper function to execute code and capture print output
def execute_and_capture_output(code):
    old_stdout = sys.stdout
    sys.stdout = new_stdout = io.StringIO()
    try:
        exec(code)
    finally:
        sys.stdout = old_stdout
    return new_stdout.getvalue()

# Main code to be executed
code_to_run = """
# Data from the problem description
mortality_control = 10
mortality_A_general = 35
mortality_B_general = 20
mortality_C_general = 10

productivity_C_infected = 50 # e.g., for Canola
productivity_C_not_infected = 30 # e.g., for Canola

# Step 1: Determine if Fungus A is a pathogen.
# A pathogen increases mortality compared to a non-infected control.
is_A_pathogen = mortality_A_general > mortality_control
print(f"Comparing Fungus A mortality to control: {mortality_A_general}% > {mortality_control}%.")
if is_A_pathogen:
    print("Conclusion: Fungus A is a pathogen.")
else:
    print("Conclusion: Fungus A is not a pathogen.")

# Step 2: Determine if Fungus B is a pathogen.
is_B_pathogen = mortality_B_general > mortality_control
print(f"\\nComparing Fungus B mortality to control: {mortality_B_general}% > {mortality_control}%.")
if is_B_pathogen:
    print("Conclusion: Fungus B is a pathogen.")
else:
    print("Conclusion: Fungus B is not a pathogen.")

# Step 3: Determine if Fungus C is a pathogen or commensal.
# A commensal does not cause harm (i.e., does not increase mortality).
is_C_pathogen = mortality_C_general > mortality_control
print(f"\\nComparing Fungus C mortality to control: {mortality_C_general}% == {mortality_control}%.")
if is_C_pathogen:
    print("Conclusion: Fungus C is a pathogen.")
else:
    print("Conclusion: Fungus C is not a pathogen, as it does not increase mortality.")
    # Check for commensal/mutualist relationship by looking at productivity
    if productivity_C_infected > productivity_C_not_infected:
        print(f"Additionally, Fungus C can increase productivity (e.g., from {productivity_C_not_infected} to {productivity_C_infected} eggs).")
        print("This supports classifying Fungus C as a commensal.")

# Final Conclusion based on the analysis
if is_A_pathogen and is_B_pathogen and not is_C_pathogen:
    print("\\nFinal Summary: The data shows that Fungus A and B are pathogens, and Fungus C is a commensal.")
"""

# Execute the code and print its output
output = execute_and_capture_output(code_to_run)
print(output)
print("<<<I>>>")