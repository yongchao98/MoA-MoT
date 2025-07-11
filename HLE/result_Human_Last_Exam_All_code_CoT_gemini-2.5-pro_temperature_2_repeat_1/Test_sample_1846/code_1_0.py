def explain_tnbc_treatment_benefit():
    """
    This function explains which population group with Triple Negative Breast Cancer (TNBC)
    shows prolonged overall survival with PD-1 inhibitors.
    """
    
    print("Based on landmark clinical trials, the benefit of adding PD-1/PD-L1 inhibitors to chemotherapy for Triple Negative Breast Cancer (TNBC) is most significant in a specific subgroup.")
    
    print("\nKey Clinical Trials:")
    print("1. KEYNOTE-355 Trial: Investigated pembrolizumab (a PD-1 inhibitor) plus chemotherapy.")
    print("2. IMpassion130 Trial: Investigated atezolizumab (a PD-L1 inhibitor) plus chemotherapy.")
    
    print("\nKey Finding for Overall Survival (OS):")
    print("Both trials demonstrated that the statistically significant improvement in overall survival was observed in the patient population whose tumors tested positive for PD-L1.")
    
    print("\nConclusion:")
    print("While some benefits might be seen in the broader 'intention-to-treat' population for other endpoints, the specific advantage of prolonged *overall survival* is robustly associated with the 'PD-L1-positive population'.")

if __name__ == "__main__":
    explain_tnbc_treatment_benefit()