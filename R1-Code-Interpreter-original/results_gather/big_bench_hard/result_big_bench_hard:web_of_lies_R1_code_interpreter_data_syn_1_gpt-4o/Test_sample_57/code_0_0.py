# Define the truth values of each person's statement
jerry_lies = True
helene_truth = jerry_lies  # Helene says Jerry lies
kandi_truth = helene_truth  # Kandi says Helene tells the truth
jaymie_lies = not kandi_truth  # Jaymie says Kandi lies
shalonda_truth = not jaymie_lies  # Shalonda says Jaymie lies

# Print the result for Shalonda
print(shalonda_truth)