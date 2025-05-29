# Define sets for Jessica's and James' circles
jessica_males = 5
jessica_females = 5
james_males = 4
james_females = 1

# We know Matilda is in both circles and is female
# So in Jessica's circle:
# - Matilda is one of the 5 females
# In James' circle:
# - Matilda is the 1 female

# Calculate Matilda's female colleagues:
# From Jessica's circle: all females except herself
matilda_females_from_jessica = jessica_females - 1

# From James' circle: no additional females (as Matilda is the only female)
matilda_females_from_james = 0

# Total female colleagues of Matilda
total_female_colleagues = matilda_females_from_jessica + matilda_females_from_james

print(f"Matilda has {total_female_colleagues} female colleagues")