# This script identifies and prints the correct statements about multichannel quantum scattering.

# Analysis Summary:
# Statement 1: False. A trivially coupled potential can lead to a non-diagonal S-matrix.
# Statement 2: True. A diagonal S-matrix implies no channel mixing, hence a diagonal potential.
# Statement 3: True. Contrapositive: a diagonal F implies a diagonal S, which implies a diagonal V (not NTC).
# Statement 4: True. Contrapositive: a diagonal S implies a diagonal V, which implies a diagonal F (not NTC).
# Statement 5: False. This contradicts statement 3.

correct_statements = [2, 3, 4]

print("The correct statements are:")
for statement_number in correct_statements:
    # This loop prints each number of the correct statements as requested.
    # The final equation is the set of correct statement numbers.
    print(statement_number)
