# Time estimates in hours
# Option A: EfficientNet with 5 species
a_train_scrape_hours = 36
a_deploy_hours = 13.8

# Option B: EfficientNet with 500 species
b_train_scrape_hours = 126
b_deploy_hours = 13.8

# Option C: ResNet with 500 species
c_train_scrape_hours = 128
c_deploy_hours = 11.8

# Option D: Manual collection
d_deploy_hours = 410

# Calculate total time for each option
total_a = a_train_scrape_hours + a_deploy_hours
total_b = b_train_scrape_hours + b_deploy_hours
total_c = c_train_scrape_hours + c_deploy_hours
total_d = d_deploy_hours

# Print the results
print("Calculating the total time for each method:")
print(f"A: {a_train_scrape_hours} hours + {a_deploy_hours} hours = {total_a} hours")
print(f"B: {b_train_scrape_hours} hours + {b_deploy_hours} hours = {total_b} hours")
print(f"C: {c_train_scrape_hours} hours + {c_deploy_hours} hours = {total_c} hours")
print(f"D: Total time = {total_d} hours")
print("\n--- Analysis ---")
print(f"Manual data collection (D) is the most time-consuming at {total_d} hours.")
print(f"Training a model for only 5 species (A) is fast at {total_a} hours, but it fails to meet the goal of identifying all pollinators.")
print(f"Training a large-scale model (B and C) takes {total_b} hours, which is the most time-efficient method that addresses the project's scope.")
print("Since methods B and C have identical total times, they represent an equally effective and time-efficient strategy compared to the others.")
