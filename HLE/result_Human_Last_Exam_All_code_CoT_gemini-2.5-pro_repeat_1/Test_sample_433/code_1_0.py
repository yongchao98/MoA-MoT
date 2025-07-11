# Define the parameters of the scenario
container_type = "single-dose container ampule"
environment = "sterile environment (e.g., ISO Class 5)"

# Explain the reasoning based on pharmaceutical standards (USP <797>)
print(f"Scenario: Determining the Beyond-Use Date (BUD) for a {container_type}.")
print(f"Opened in a: {environment}.")
print("---------------------------------------------------------------------")
print("Guideline based on United States Pharmacopeia (USP) standards:")
print("1. An ampule is a sealed glass container that is broken open to access its contents.")
print("2. Once an ampule is opened, it cannot be resealed to protect it from environmental contamination.")
print("3. Because it cannot be resealed, sterility cannot be maintained after opening.")

# State the final conclusion for the BUD
print("\nConclusion:")
print("The contents of the ampule must be used immediately after opening.")
print("Any unused portion must be discarded and cannot be stored for any length of time.")
print("\nTherefore, the official Beyond-Use Date (BUD) is: Immediate Use.")
