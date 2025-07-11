# Step 1: Define the numerical values for the letters in "Кома" based on the Russian alphabet.
# К is the 12th letter, О is the 16th, М is the 14th, and А is the 1st.
k_val = 12
o_val = 16
m_val = 14
a_val = 1

print("Decoding 'Кома' based on Russian alphabet positions:")
print(f"К = {k_val}")
print(f"О = {o_val}")
print(f"М = {m_val}")
print(f"А = {a_val}")
print("-" * 20)

# Step 2: Calculate the Latitude using the derived formula.
# Latitude = К + О + М + А + О
latitude = k_val + o_val + m_val + a_val + o_val
print("Calculating Latitude:")
print(f"Equation: {k_val} + {o_val} + {m_val} + {a_val} + {o_val} = {latitude}")
print("-" * 20)

# Step 3: Calculate the Longitude using the derived formula.
# Longitude = М * (О - К)
longitude = m_val * (o_val - k_val)
print("Calculating Longitude:")
print(f"Equation: {m_val} * ({o_val} - {k_val}) = {longitude}")
print("-" * 20)

# Step 4: Output the final coordinates and identify the location.
print(f"The final coordinates are approximately {latitude}° N, {longitude}° E.")
print("These coordinates point to Perm Krai.")
