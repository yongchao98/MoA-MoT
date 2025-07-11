# Time estimates for each method in hours
time_a = {'train': 36, 'deploy': 13.8}
time_b = {'train': 126, 'deploy': 13.8}
time_c = {'train': 128, 'deploy': 11.8}
time_d = {'deploy': 410}

# The goal is to identify ALL pollinators, so Option A (5 species) is insufficient.

# Calculate total time for Option B
total_time_b = time_b['train'] + time_b['deploy']
print(f"Total time for Option B (EfficientNet with 500 species):")
print(f"{time_b['train']} hours (training) + {time_b['deploy']} hours (deployment) = {total_time_b} hours")
print("-" * 30)

# Calculate total time for Option C
total_time_c = time_c['train'] + time_c['deploy']
print(f"Total time for Option C (ResNet with 500 species):")
print(f"{time_c['train']} hours (training) + {time_c['deploy']} hours (deployment) = {total_time_c} hours")
print("-" * 30)

# Total time for Option D
total_time_d = time_d['deploy']
print(f"Total time for Option D (Manual collection):")
print(f"{total_time_d} hours")
print("-" * 30)

# Conclusion
print("\nComparison:")
print(f"Option B takes {total_time_b} hours.")
print(f"Option C takes {total_time_c} hours.")
print(f"Option D takes {total_time_d} hours.")
print("\nBoth Option B and Option C require significantly less time than Option D.")
print("Since both B and C are the fastest and require the same amount of total time, they are the easiest methods.")
print("Therefore, the best answer choice is the one that includes both B and C.")
