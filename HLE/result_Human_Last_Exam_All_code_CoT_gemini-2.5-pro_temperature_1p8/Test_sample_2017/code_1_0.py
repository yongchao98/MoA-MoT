# W1: The non-classical word in Passage 1 is εἶμαι (eimai).
# C1: The correct classical Attic form is εἰμί (eimi).
# P1: The form εἶμαι is found in Koine Greek and is the standard in Demotic (Modern) Greek.
W1 = "εἶμαι"
C1 = "εἰμί"
P1 = "KoineDemotic"

# W2: The non-classical word in Passage 2 is ὄρνιθα (ornitha).
# C2: The correct classical Attic form for the accusative singular is ὄρνιν (ornin).
# P2: The form ὄρνιθα is found in Homeric Greek and becomes standard in Koine Greek.
W2 = "ὄρνιθα"
C2 = "ὄρνιν"
P2 = "HomericKoine"

# Combine the results into the specified comma-separated format.
final_string = f"{W1},{C1},{P1},{W2},{C2},{P2}"

# Print the final string.
print(final_string)