# Parameters from the problem description
c = 0.95  # Consistency level of decentralized identifiers
b = 3     # Branching factor of semantic version control

# --- Model anad Calculation ---

# 1. Model the theoretical maximum value for each FAIR component
# Findability and Accessibility are primarily limited by identifier consistency.
f_max = c
a_max = c

# Interoperability and Reusability are limited by identifier consistency
# and also degraded by the complexity of semantic branching.
# We model this degradation as a factor of 1/b.
i_max = c / b
r_max = c / b

# 2. Calculate the overall FAIR score R as the average of the components.
R_max = (f_max + a_max + i_max + r_max) / 4


# --- Output the results ---

print("Step 1: Define the model for the FAIR score R.")
print("R is the arithmetic mean of its four components:")
print("R = (f + a + i + r) / 4")
print("-" * 30)

print("Step 2: Calculate the maximum value for each component based on given parameters.")
print(f"Identifier consistency 'c' = {c}")
print(f"Semantic branching factor 'b' = {b}\n")

print(f"Findability (f_max) = c = {f_max}")
print(f"Accessibility (a_max) = c = {a_max}")
print(f"Interoperability (i_max) = c / b = {c} / {b} = {i_max}")
print(f"Reusability (r_max) = c / b = {c} / {b} = {r_max}")
print("-" * 30)

print("Step 3: Calculate the theoretical maximum FAIR score R_max.")
print("The final equation is:")
print(f"R_max = (f_max + a_max + i_max + r_max) / 4")
print("\nPlugging in the numbers:")
# The f-string below shows each number in the final equation.
print(f"R_max = ({f_max} + {a_max} + {i_max} + {r_max}) / 4")

# Display the final calculated value
print(f"\nThe theoretical maximum value of R is: {R_max}")

# The final answer in the requested format
print(f"\n<<<{R_max}>>>")