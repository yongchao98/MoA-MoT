# Step 1: Identify the non-classical word, its classical version, and its period from Passage 1.
# The word "εἶμαι" is the form of "I am" used in Koine and Demotic Greek. The classical form is "εἰμί".
W1 = "εἶμαι"
C1 = "εἰμί"
P1 = "KoineDemotic"

# Step 2: Identify the non-classical word, its classical version, and its period from Passage 2.
# The word "ὄρνιθα" is the Homeric accusative singular of "ὄρνις" (bird). The classical Attic form is "ὄρνιν".
W2 = "ὄρνιθα"
C2 = "ὄρνιν"
P2 = "Homeric"

# Step 3: Format the results into a single comma-separated string as specified.
# The required format is "$W_1$,$C_1$,$P_1$,$W_2$,$C_2$,$P_2$".
final_result = f"{W1},{C1},{P1},{W2},{C2},{P2}"

# Step 4: Print the final formatted string.
print(final_result)