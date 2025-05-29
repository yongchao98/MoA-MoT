# Jessica's circle (excluding Matilda)
jessica_males = 6
jessica_females = 1  # This includes Matilda

# Michael's circle (excluding Matilda)
michael_males = 5
michael_females = 6  # This includes Matilda

# Calculate Matilda's female colleagues
# From Jessica's circle:
matilda_females = 1  # Jessica herself

if jessica_females > 0:
    matilda_females += (jessica_females - 1)  # The other female (if any) minus Matilda

# From Michael's circle:
matilda_females += (michael_females - 1)  # All females minus Matilda herself

print(f"Matilda has {matilda_females} female colleagues")