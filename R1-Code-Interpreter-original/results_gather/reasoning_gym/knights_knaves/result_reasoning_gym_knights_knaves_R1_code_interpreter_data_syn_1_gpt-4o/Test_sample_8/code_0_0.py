# Define the truth values for Victoria and Grace
victoria_is_angel = True
grace_is_angel = False

# Check Victoria's statement: "Victoria is an angel or Grace is an angel"
victoria_statement = victoria_is_angel or grace_is_angel

# Check Grace's statement: "Grace is an angel if and only if Victoria is a devil"
grace_statement = (grace_is_angel == (not victoria_is_angel))

# Print the results
print("Victoria's statement is", victoria_statement)
print("Grace's statement is", grace_statement)