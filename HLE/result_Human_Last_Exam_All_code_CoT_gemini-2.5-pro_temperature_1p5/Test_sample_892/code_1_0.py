def medical_disclaimer():
  """
  This function prints a disclaimer about the limitations of AI in providing medical advice.
  """
  patient_age = 57
  creatinine_level = 2.1

  print("Disclaimer: I am an AI assistant and not a medical professional.")
  print("The following information is for illustrative purposes and does not constitute a medical diagnosis.")
  print("-" * 50)
  print(f"The clinical case presented involves a {patient_age}-year-old woman with multiple significant findings:")
  print("- A history of COPD")
  print("- Symptoms of dyspnea, chronic cough, and acid reflux")
  print("- A vertebral mass on CT scan")
  print(f"- Elevated blood creatinine of {creatinine_level}")
  print("-" * 50)
  print("These findings require careful evaluation by a qualified healthcare provider.")
  print("Making a diagnosis is a complex process that cannot be done by an AI.")
  print("\nPlease consult a doctor for any medical concerns.")

medical_disclaimer()