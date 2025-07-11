import textwrap

def explain_and_demonstrate_thooft_anomaly_matching():
    """
    This script explains and demonstrates the physical implication of the
    't Hooft anomaly matching condition.
    """
    
    # Step 1: Define a high-energy (UV) theory.
    # We imagine a theory like QCD with fundamental quarks that have a global chiral symmetry.
    uv_theory_name = "High-Energy QCD-like Theory (UV)"
    
    # Step 2: This theory has a 't Hooft anomaly.
    # An anomaly is a quantum mechanical breaking of a classical symmetry.
    # The anomaly's value is a robust, calculable number that doesn't change
    # as we look at the theory at different energy scales.
    # Let's say we calculate it from the fundamental quarks and get a specific value.
    uv_anomaly_value = 1.0  # This value is fixed and must be preserved.

    print(textwrap.dedent(f"""
    ===============================================================
    't Hooft Anomaly Matching: A Demonstration
    ===============================================================
    
    1. The High-Energy (UV) Theory
    ---------------------------------
    We start with a theory: '{uv_theory_name}'
    This theory has fundamental particles (like quarks).
    It possesses a global symmetry which has a 't Hooft Anomaly.
    
    The calculated anomaly value in the UV is: {uv_anomaly_value}
    """))

    # Step 3: State the 't Hooft Anomaly Matching Condition.
    # The anomaly calculated in the low-energy (IR) effective theory must
    # exactly match the one from the UV theory.
    # The Final Equation is: Anomaly(UV) == Anomaly(IR)
    
    print(textwrap.dedent(f"""
    2. The 't Hooft Matching Condition
    ---------------------------------
    The rule states that any valid low-energy (IR) theory that describes
    the same physics must reproduce this exact anomaly.
    
    The Governing Equation: Anomaly(UV) == Anomaly(IR)
    Our Target Value: {uv_anomaly_value}
    """))
    
    # Step 4: Propose several possible low-energy (IR) theories.
    # These theories describe physics in terms of composite particles (like hadrons).
    # We must check if they satisfy the condition.
    
    possible_ir_theories = [
        {
            "name": "IR Theory #1: Trivial Gapped Theory",
            "description": "A theory where all composite particles are heavy and the symmetry is preserved.",
            "anomaly_value": 0.0, # Heavy particles don't contribute to the anomaly.
        },
        {
            "name": "IR Theory #2: Chiral Symmetry Breaking",
            "description": "A theory with massless Goldstone bosons (pions) from broken symmetry.",
            "anomaly_value": 1.0, # The pions' interactions correctly reproduce the anomaly.
        },
        {
            "name": "IR Theory #3: Incorrect Massless Fermions",
            "description": "A theory with massless composite fermions with the wrong quantum numbers.",
            "anomaly_value": -0.5, # These fermions generate a different, incorrect anomaly.
        }
    ]
    
    print(textwrap.dedent("""
    3. Testing Proposed Low-Energy (IR) Theories
    ----------------------------------------------
    We now test which of our proposed IR theories are valid."""))

    # Step 5: Test each proposed IR theory.
    for theory in possible_ir_theories:
        print(f"\n--- Checking: {theory['name']} ---")
        print(f"Description: {theory['description']}")
        
        # This is the check: does the IR anomaly match the UV anomaly?
        ir_anomaly = theory['anomaly_value']
        
        print(f"Checking Equation: {uv_anomaly_value} == {ir_anomaly}")
        
        if uv_anomaly_value == ir_anomaly:
            print(">>> RESULT: VALID. The anomaly matches. This theory is a possible description of the low-energy physics.")
        else:
            print(">>> RESULT: INVALID. The anomaly does not match. This theory is ruled out by the 't Hooft condition.")
            
    # Step 6: Conclusion
    print(textwrap.dedent("""
    ===============================================================
    4. Conclusion: The Physical Implication
    ===============================================================
    As demonstrated, the anomaly matching condition is not just a mathematical curiosity.
    It acts as a powerful **CONSTRAINT** that any proposed low-energy effective theory
    must satisfy. It connects the deep UV physics to the observable IR physics in a
    non-trivial, non-perturbative way, allowing us to rule out entire classes of
    theories that would otherwise seem plausible. This is its key implication.
    """))

if __name__ == '__main__':
    explain_and_demonstrate_thooft_anomaly_matching()