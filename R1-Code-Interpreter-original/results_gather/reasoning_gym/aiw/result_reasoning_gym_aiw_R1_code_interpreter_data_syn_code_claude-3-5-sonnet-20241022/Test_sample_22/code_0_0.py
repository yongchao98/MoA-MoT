# Let's create sets for Jennifer's and Joseph's circles
# We know Matilda is in both circles

# Jennifer's circle (excluding Jennifer herself)
jennifer_males = 5
jennifer_females = 3  # including Matilda

# Joseph's circle (excluding Joseph himself)
joseph_males = 4
joseph_females = 4  # including Matilda

# Since Matilda is in both circles and is female
# Let's calculate Matilda's female colleagues:

# In Jennifer's circle:
matilda_females_from_jennifer = jennifer_females - 1  # excluding Matilda herself

# In Joseph's circle:
matilda_females_from_joseph = joseph_females - 1  # excluding Matilda herself

# Total unique female colleagues of Matilda
total_matilda_females = matilda_females_from_jennifer + matilda_females_from_joseph

print(f"Matilda has {total_matilda_females} female colleagues")