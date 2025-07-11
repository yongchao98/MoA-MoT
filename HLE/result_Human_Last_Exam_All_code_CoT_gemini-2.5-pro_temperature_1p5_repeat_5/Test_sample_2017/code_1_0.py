# Word from passage 1 that is not classical
W1 = "εἶμαι"
# Correct classical version for W1
C1 = "εἰμί"
# Period(s) to which W1 belongs
P1 = "KoineByzantineDemotic"

# Word from passage 2 that is not classical
W2 = "ὄρνιθα"
# Correct classical version for W2
C2 = "ὄρνιν"
# Period(s) to which W2 belongs
P2 = "HomericKoine"

# Combine all the parts into the final string, separated by commas.
final_answer = ",".join([W1, C1, P1, W2, C2, P2])

# Print the final result.
print(final_answer)