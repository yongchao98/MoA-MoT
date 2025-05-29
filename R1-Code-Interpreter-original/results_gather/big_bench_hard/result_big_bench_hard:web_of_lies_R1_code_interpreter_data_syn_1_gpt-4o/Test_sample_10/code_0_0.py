# Define the truth values
jerry_truth = False  # Jerry lies
millicent_truth = jerry_truth  # Millicent says Jerry tells the truth
sal_truth = millicent_truth  # Sal says Millicent tells the truth
amberly_truth = sal_truth  # Amberly says Sal tells the truth
shenna_truth = not amberly_truth  # Shenna says Amberly lies

# Print the result for Shenna
print(shenna_truth)