# This script analyzes the entities involved in the real estate dispute.

# Number of primary companies involved in the dispute
num_companies = 2 # RE1 and RE2

# Number of distinct residential cities of the principals
# RE1 principals: Richmond Hill, King City, Markham
# RE2 principals: Mississauga, Vaughan
num_principal_locations = 5 

# Number of commercial real estate properties mentioned
# Ottawa, North Bay, Sudbury, Markham, Niagara Falls, Oshawa
num_properties = 6

# Total count of these specific entities
total_entities = num_companies + num_principal_locations + num_properties

print("This script quantifies the key parties and properties in the dispute.")
print("The analysis involves a simple summation of these entities.")
print("Equation based on the entities described:")
print(num_companies, "+", num_principal_locations, "+", num_properties, "=", total_entities)
