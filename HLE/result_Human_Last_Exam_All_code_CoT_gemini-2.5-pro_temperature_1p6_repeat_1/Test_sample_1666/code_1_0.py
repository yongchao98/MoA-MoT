import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
string_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = string_buffer

print("Here is a step-by-step legal analysis of the four questions:")

# Question 1: Bryan's Restrictive Covenants (CEO)
print("\n1) Bryan's Employment Agreement (CEO):")
print("   - Context: The agreement is ancillary to the sale of a business, a scenario where courts are more likely to enforce restrictive covenants to protect the buyer's purchased goodwill.")
print("   - Ontario Law: The statutory ban on non-competes does not apply to executives (like a CEO) or in the context of a business sale. Bryan meets both exceptions.")
print("   - Conclusion: The 1-year non-solicitation and 6-month non-competition clauses are likely reasonable and enforceable.")

# Question 2: Ryan's Restrictive Covenants (Shift Manager)
print("\n2) Ryan's Employment Agreement (Shift Manager):")
print("   - Context: While also part of the sale, the reasonableness of the clauses is judged against his specific role.")
print("   - Analysis: For a non-executive Shift Manager, a non-competition clause covering all of Ontario is almost certainly an unreasonable restraint of trade. A non-solicitation clause is also difficult to enforce for a role with minimal client contact.")
print("   - Conclusion: Both the non-competition and non-solicitation clauses are likely unenforceable.")

# Question 3: New Employees' Employment Agreements
print("\n3) New Employees' Employment Agreements:")
print("   - Context: These are new hires in manufacturing roles, not executives or part of the sale.")
print("   - Ontario Law: The statutory ban on non-competition clauses applies directly, making those clauses void.")
print("   - Severability: The presence of a void clause does not invalidate the entire agreement. A court will sever the bad clause and enforce the rest of the contract.")
print("   - Conclusion: The agreements are valid, but the non-competition clauses within them are unenforceable.")

# Question 4: The Pickup Truck
print("\n4) The Promise to Transfer the Pickup Truck:")
print("   - Context: Bryan's promise was made 'out of appreciation' and was not part of any formal agreement.")
print("   - Legal Principle: An enforceable contract requires consideration (an exchange of value). Stan provided no consideration for the promise of the truck.")
print("   - Conclusion: This was a gratuitous promise (a promise to make a gift) and is not legally enforceable. Bryan is not required to transfer the truck.")

# Final Synthesis and Answer
print("\n---")
print("Summary of Analysis:")
print(" - Bryan’s clauses are VALID and ENFORCEABLE.")
print(" - Ryan’s clauses are NOT ENFORCEABLE.")
print(" - The new employees' agreements are VALID, but their non-competition clauses are UNENFORCEABLE.")
print(" - Bryan is NOT REQUIRED to transfer the truck.")
print("\nThis combination of outcomes corresponds directly to Answer Choice A.")
print("\n<<<A>>>")

# Get the content from the buffer
output = string_buffer.getvalue()
# Restore original stdout
sys.stdout = original_stdout
# Print the captured output
print(output)