import textwrap

def analyze_newborn_case():
    """
    Analyzes the clinical vignette and identifies the most likely cause.
    """

    # Clinical findings from the problem description
    weight_lb = 12
    weight_oz = 1
    oxygen_saturation = 89
    
    # Print the analysis step-by-step
    print("Analyzing the Clinical Vignette:")
    print("-" * 35)

    # Step 1: Analyze Macrosomia
    analysis_macrosomia = f"""
    1.  Macrosomia (Large for Gestational Age): The newborn's weight is {weight_lb}-lb {weight_oz}-oz.
        This is significantly large for a newborn. The most common cause of fetal macrosomia is poorly controlled maternal diabetes.
    """
    print(textwrap.dedent(analysis_macrosomia))

    # Step 2: Analyze Respiratory Distress
    analysis_resp_distress = f"""
    2.  Respiratory Distress: The decreased oxygen saturation of {oxygen_saturation}% and a fluid-filled lung density indicate respiratory distress.
        Infants of diabetic mothers are at increased risk for Neonatal Respiratory Distress Syndrome (NRDS) because high fetal insulin levels can impair surfactant production.
    """
    print(textwrap.dedent(analysis_resp_distress))

    # Step 3: Analyze Congenital Anomaly
    analysis_anomaly = """
    3.  Congenital Anomaly: The presence of micrognathia (small jaw) is a congenital defect.
        Maternal diabetes is a known teratogen, increasing the risk of various anomalies, including craniofacial defects.
    """
    print(textwrap.dedent(analysis_anomaly))

    # Final Conclusion
    conclusion = """
    Conclusion:
    The combination of all three key findings (macrosomia, respiratory distress, and a congenital anomaly) points strongly to maternal diabetes as the single unifying underlying cause.
    """
    print(conclusion)

if __name__ == "__main__":
    analyze_newborn_case()