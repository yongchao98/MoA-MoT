# In Passage 1, the non-classical word is εἶμαι (I am). The correct classical form is εἰμί.
# The form εἶμαι is found in Koine, Byzantine, and Modern (Demotic) Greek.
w1 = "εἶμαι"
c1 = "εἰμί"
p1 = "KoineByzantineDemotic"

# In Passage 2, the non-classical word is ὄρνιθα (bird, accusative singular).
# The correct classical Attic prose form is ὄρνιν. The form ὄρνιθα is characteristic of Homeric Greek.
w2 = "ὄρνιθα"
c2 = "ὄρνιν"
p2 = "Homeric"

# We format the final answer as a comma-separated string.
final_answer = f"{w1},{c1},{p1},{w2},{c2},{p2}"

# We print the final result.
print(final_answer)
