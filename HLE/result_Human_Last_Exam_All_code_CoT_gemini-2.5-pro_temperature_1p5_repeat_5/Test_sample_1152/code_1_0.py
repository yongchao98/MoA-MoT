# Data extracted from the problem description
# All mortality rates are in percent (%)
control_mortality = 10
fungus_a_max_mortality = 35
fungus_b_mortality = 20
fungus_c_mortality = 10

# Step 1: Analyze if Fungus A is a pathogen
print("--- Analysis of Fungus A ---")
print(f"The mortality rate for non-infected honeybees (control) is {control_mortality}%.")
print(f"The highest mortality rate observed for bees infected with Fungus A is {fungus_a_max_mortality}%.")
# A pathogen increases the host's mortality rate. Let's check the condition.
print(f"To check if Fungus A is a pathogen, we compare its mortality rate to the control: {fungus_a_max_mortality} > {control_mortality}")
if fungus_a_max_mortality > control_mortality:
    print("Result: True. Fungus A significantly increases mortality, so it is a pathogen.\n")
else:
    print("Result: False. Fungus A does not increase mortality, so it is not a pathogen.\n")

# Step 2: Analyze if Fungus B is a pathogen
print("--- Analysis of Fungus B ---")
print(f"The mortality rate for non-infected honeybees (control) is {control_mortality}%.")
print(f"The mortality rate observed for bees infected with Fungus B is {fungus_b_mortality}%.")
# Check if Fungus B increases the host's mortality rate.
print(f"To check if Fungus B is a pathogen, we compare its mortality rate to the control: {fungus_b_mortality} > {control_mortality}")
if fungus_b_mortality > control_mortality:
    print("Result: True. Fungus B increases mortality, so it is a pathogen.\n")
else:
    print("Result: False. Fungus B does not increase mortality, so it is not a pathogen.\n")

# Step 3: Analyze if Fungus C is a pathogen or a commensal
print("--- Analysis of Fungus C ---")
print(f"The mortality rate for non-infected honeybees (control) is {control_mortality}%.")
print(f"The mortality rate observed for bees infected with Fungus C is {fungus_c_mortality}%.")
# A commensal organism does not harm the host, meaning mortality does not increase.
print(f"To check if Fungus C is a commensal, we see if its mortality rate is equal to the control: {fungus_c_mortality} == {control_mortality}")
if fungus_c_mortality <= control_mortality:
    print("Result: True. Fungus C does not increase mortality, so it is a commensal.\n")
else:
    print("Result: False. Fungus C increases mortality, so it is a pathogen.\n")

# Step 4: Final Conclusion
print("--- Final Conclusion ---")
print("Our analysis shows that Fungus A and Fungus B are pathogens because they cause mortality rates of 35% and 20% respectively, both higher than the 10% control rate.")
print("Fungus C is a commensal because its associated mortality rate of 10% is the same as the control rate, indicating it does not harm the bees.")
print("The statement 'Fungus A and B are pathogens. Fungus C is a commensal.' is fully supported by the data.")