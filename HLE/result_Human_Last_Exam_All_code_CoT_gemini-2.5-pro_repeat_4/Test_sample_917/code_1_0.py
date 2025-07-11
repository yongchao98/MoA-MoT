# Data for each method: {Name: [Training_Time, Deployment_Time]}
methods = {
    "A. EfficientNet (5 species)": [36, 13.8],
    "B. EfficientNet (500 species)": [126, 13.8],
    "C. ResNet (500 species)": [128, 11.8],
    "D. Manual": [0, 410] # Manual method has no "training" time, only "deployment"
}

print("Calculating total time for each method:\n")

for method_name, times in methods.items():
    training_time = times[0]
    deployment_time = times[1]
    total_time = training_time + deployment_time
    
    # Print the equation and the result
    if training_time > 0:
        print(f"Method {method_name}:")
        print(f"{training_time} + {deployment_time} = {total_time:.1f} hours")
    else:
        # For manual method
        print(f"Method {method_name}:")
        print(f"{deployment_time} = {deployment_time:.1f} hours")
    
    # Add a note about feasibility
    if "5 species" in method_name:
        print("Note: This method is fast but fails to meet the goal of identifying all pollinators.\n")
    elif "500 species" in method_name:
        print("Note: This method is feasible and significantly faster than manual collection.\n")
    else:
        print("Note: This is the most time-consuming method.\n")


print("Conclusion:")
print("Method A is not suitable as it won't identify most of the pollinators.")
print("Method D is too slow (410 hours).")
print("Methods B and C are both the most efficient and feasible options, taking the same total time (139.8 hours).")
print("Therefore, both B and C represent the easiest valid methods.")
