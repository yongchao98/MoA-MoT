import sys
from io import StringIO

class EmployerAnalyzer:
    """Analyzes employer compliance with Ontario employment laws as of Jan 2, 2023."""

    EMPLOYEE_THRESHOLD = 25

    def __init__(self, name, description, emp_count_jan1_2022, has_disconnect_policy, has_monitoring_policy, distributes_policies):
        self.name = name
        self.description = description
        self.emp_count_jan1_2022 = emp_count_jan1_2022
        self.has_disconnect_policy = has_disconnect_policy
        self.has_monitoring_policy = has_monitoring_policy
        self.distributes_policies = distributes_policies
        self.is_compliant = True
        self.reason = []

    def check_compliance(self):
        """Checks compliance for disconnecting from work and electronic monitoring policies."""
        print(f"--- Analyzing Employer {self.name} ---")

        # The employee count is determined on Jan 1, 2022, for policies required in 2022/2023.
        requires_policy = self.emp_count_jan1_2022 >= self.EMPLOYEE_THRESHOLD
        
        print(f"Employee count on Jan 1, 2022: {self.emp_count_jan1_2022}")
        print(f"Employee threshold for requiring a policy: {self.EMPLOYEE_THRESHOLD}")

        if not requires_policy:
            print(f"Result: As {self.emp_count_jan1_2022} is less than {self.EMPLOYEE_THRESHOLD}, this employer was NOT required to have these policies.")
            self.reason.append(f"Compliant because they are below the {self.EMPLOYEE_THRESHOLD} employee threshold.")
        else:
            print(f"Result: As {self.emp_count_jan1_2022} is greater than or equal to {self.EMPLOYEE_THRESHOLD}, this employer was REQUIRED to have both policies.")
            # Check for Disconnecting From Work Policy
            if not self.has_disconnect_policy:
                self.is_compliant = False
                self.reason.append("NON-COMPLIANT: Failed to develop a required 'Disconnecting from Work' policy.")
            
            # Check for Electronic Monitoring Policy
            if not self.has_monitoring_policy:
                self.is_compliant = False
                self.reason.append("NON-COMPLIANT: Failed to develop a required 'Electronic Monitoring' policy.")

            # Check for Distribution of required policies
            if not self.distributes_policies:
                self.is_compliant = False
                self.reason.append("NON-COMPLIANT: Failed to distribute required policies to employees.")
        
        if not self.reason:
            self.reason.append("Compliant as all legal requirements are met.")

        status = "COMPLIANT" if self.is_compliant else "NOT IN COMPLIANCE"
        print(f"Final Status: {status}\n")
        return self.is_compliant, self.name, self.reason[0]

def main():
    """Main function to analyze all employers."""
    
    # Note on Employer B: The Jan 1, 2022 count is not given. We assume they were below the threshold,
    # as the current count is 23. If they were above the threshold, their refusal to distribute policies
    # would also be a violation. However, Employer D presents a clearer, fact-based violation.
    employers_data = [
        ("A", "Has 20 employees, no disconnect policy.", 20, False, True, True),
        ("B", "Has 23 employees, refuses to distribute policies.", 23, True, True, False), # Assuming 23 was the count on Jan 1, 2022 for this analysis
        ("C", "Has 1,000 employees, has and distributes both policies.", 1000, True, True, True),
        ("D", "Has 30 employees, has disconnect policy but no monitoring policy.", 30, True, False, True),
        ("E", "Has 22 employees, no disconnect policy.", 22, False, True, True)
    ]
    
    non_compliant_employer = None
    
    # Backup the original stdout
    original_stdout = sys.stdout
    # Create a string buffer to capture the output
    captured_output = StringIO()
    # Redirect stdout to the buffer
    sys.stdout = captured_output

    for data in employers_data:
        analyzer = EmployerAnalyzer(*data)
        is_compliant, name, reason = analyzer.check_compliance()
        if not is_compliant and non_compliant_employer is None:
            non_compliant_employer = {"name": name, "reason": reason}

    # Restore the original stdout
    sys.stdout = original_stdout
    
    # Get the captured output as a string
    output_string = captured_output.getvalue()
    print(output_string)

    if non_compliant_employer:
        print("-----------------------------------------")
        print(f"The employer not in compliance is: Employer {non_compliant_employer['name']}")
        print(f"Reasoning: {non_compliant_employer['reason']}")
        print("-----------------------------------------")
    else:
        print("All employers appear to be in compliance based on the data.")

if __name__ == "__main__":
    main()
