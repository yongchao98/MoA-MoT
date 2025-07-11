# In Passage 1, the verb for "I am" is given as εἶμαι.
# The standard classical form is εἰμί. The form εἶμαι is the standard form
# in Koine Greek and subsequent periods.
W1 = "εἶμαι"
C1 = "εἰμί"
P1 = "KoineByzantineKatharevousaDemotic"

# In Passage 2, the accusative singular of "bird" is given as ὄρνιθα.
# The classical third-declension noun ὄρνις (gen. ὄρνιθος) has the accusative ὄρνιν.
# The form ὄρνιθα comes from a later re-declension of the noun as a first-declension noun,
# a change that began in the Koine period.
W2 = "ὄρνιθα"
C2 = "ὄρνιν"
P2 = "KoineByzantineDemotic"

# Combine the results into the specified format.
result = f"{W1},{C1},{P1},{W2},{C2},{P2}"

print(result)