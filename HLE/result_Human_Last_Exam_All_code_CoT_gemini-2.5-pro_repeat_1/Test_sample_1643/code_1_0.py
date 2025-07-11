def diagnose_aps():
    """
    Analyzes patient data to determine if she meets the criteria for
    Antiphospholipid Syndrome (APS) based on the revised Sapporo (Sydney) criteria.
    """

    # Patient Data
    patient_info = {
        "age": 34,
        "bmi": 32,
        "vte_events": [
            {"type": "Calf DVT", "age": 18, "trigger": "Plane travel"},
            {"type": "Intermediate-risk PE", "age": 25, "trigger": "Pregnancy"},
            {"type": "PE", "age": 34, "trigger": "Spontaneous"}
        ],
        "labs_3_months_ago": {
            "anti_b2gp1_igm": 41,
            "anticardiolipin_igm": 32,
            "ptt_la_ratio": 1.19,
            "drvv_ratio": 1.44
        },
        "labs_today": {
            "anti_b2gp1_igm": 29,
            "anti_b2gp1_igg": 21,
            "anticardiolipin_igm": 47,
            "ptt_la_ratio": 1.17,
            "drvv_ratio": 1.51
        }
    }

    # --- Diagnostic Steps ---
    print("Step 1: Assessing Clinical Criteria for APS")
    print("---------------------------------------------")
    print("The APS diagnostic criteria require at least one clinical criterion.")
    print("Clinical criteria include vascular thrombosis (venous or arterial).")
    print(f"This patient has a history of 3 venous thromboembolic (VTE) events:")
    print(f"  1. Calf DVT at age 18")
    print(f"  2. Pulmonary Embolism (PE) at age 25")
    print(f"  3. Another PE 4 months ago")
    print("Conclusion: The patient clearly meets the clinical criterion for vascular thrombosis.\n")

    print("Step 2: Assessing Laboratory Criteria for APS")
    print("---------------------------------------------")
    print("The criteria require at least one laboratory marker to be positive on two occasions, at least 12 weeks apart.")
    print("The lab tests were performed 3 months ago and today, satisfying the 12-week interval.\n")
    print("We assess three types of antiphospholipid antibodies:")

    # Lupus Anticoagulant (LA) assessment
    print("1. Lupus Anticoagulant (LA):")
    print(f"   - 3 months ago: The dRVVT ratio was {patient_info['labs_3_months_ago']['drvv_ratio']} (Normal < 1.2), which is positive.")
    print(f"   - Today: The dRVVT ratio is {patient_info['labs_today']['drvv_ratio']} (Normal < 1.2), which is persistently positive.")
    print("   Conclusion: The Lupus Anticoagulant criterion is met based on the persistent positive dRVVT.\n")

    # Anticardiolipin (aCL) and Anti-β2GP1 assessment
    print("2. Anticardiolipin (aCL) and Anti-β2-Glycoprotein-I (aβ2GPI) Antibodies:")
    print(f"   - Anticardiolipin IgM: Positive 3 months ago ({patient_info['labs_3_months_ago']['anticardiolipin_igm']} UI/L, N < 20) and positive today ({patient_info['labs_today']['anticardiolipin_igm']} UI/L, N < 20).")
    print(f"   - Anti-β2GP1 IgM: Positive 3 months ago ({patient_info['labs_3_months_ago']['anti_b2gp1_igm']} UI/L, N < 20) and positive today ({patient_info['labs_today']['anti_b2gp1_igm']} UI/L, N < 20).")
    print("   - Anti-β2GP1 IgG: Became positive today ({patient_info['labs_today']['anti_b2gp1_igg']} UI/L, N < 20).")
    print("   Conclusion: These results show persistent positivity, further supporting the diagnosis.\n")

    print("Step 3: Final Diagnosis")
    print("-------------------------")
    print("The patient meets both the clinical and laboratory criteria for Antiphospholipid Syndrome:")
    print("- Clinical Criterion: Met, due to multiple episodes of venous thrombosis (3 events).")
    print("- Laboratory Criterion: Met, due to a persistently positive Lupus Anticoagulant (dRVVT ratio of 1.44 and 1.51) and persistently positive anticardiolipin and anti-β2GP1 antibodies over a 3-month period.")
    print("\nTherefore, the patient can be diagnosed with Antiphospholipid Syndrome.")

    final_answer = "Yes"
    print(f"\nDoes this patient categorize as having antiphospholipid syndrome?")
    print(f'<<<{final_answer}>>>')

if __name__ == '__main__':
    diagnose_aps()