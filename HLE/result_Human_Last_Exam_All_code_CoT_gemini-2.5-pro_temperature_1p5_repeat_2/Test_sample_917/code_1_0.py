# Time estimates for each method in hours
# [training_time, deployment_time]
methods = {
    'A': [36, 13.8],
    'B': [126, 13.8],
    'C': [128, 11.8],
    'D': [0, 410] # Manual method has no "training", only "deployment"
}

# Calculate and print the total time for each viable method
print("Calculating total time for each viable method:")

# Method B: EfficientNet with 500 species
time_b_train = methods['B'][0]
time_b_deploy = methods['B'][1]
total_time_b = time_b_train + time_b_deploy
print(f"Method B total time = {time_b_train} + {time_b_deploy} = {total_time_b} hours")

# Method C: ResNet with 500 species
time_c_train = methods['C'][0]
time_c_deploy = methods['C'][1]
total_time_c = time_c_train + time_c_deploy
print(f"Method C total time = {time_c_train} + {time_c_deploy} = {total_time_c} hours")

# Method D: Manual collection
time_d_deploy = methods['D'][1]
total_time_d = time_d_deploy
print(f"Method D total time = {total_time_d} hours")

print("\nComparison:")
print(f"The automated methods (B and C) take {total_time_b} hours.")
print(f"The manual method (D) takes {total_time_d} hours.")

if total_time_b < total_time_d and total_time_c < total_time_d:
    print("Both automated methods B and C are significantly faster and therefore easier than the manual method D.")
    if total_time_b == total_time_c:
        print("Methods B and C are equally easy in terms of total time.")
