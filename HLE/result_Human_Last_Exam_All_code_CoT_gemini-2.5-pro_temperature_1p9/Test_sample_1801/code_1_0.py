# The task is to identify two Nobel laureates from a specific commercial.
# This is a factual knowledge question.
# I will store the names of the laureates in a list and then print them in a formatted string.

laureates = ["Paul Krugman", "Myron Scholes"]
company = "Enron"
year = 2000
tagline = "Ask Why"

print(f"The two Nobel laureates featured in the infamous {year} commercial for {company} with the tagline '{tagline}' were:")
print(f"- {laureates[0]}")
print(f"- {laureates[1]}")