import math

# Step 1: Identify the number X.
# X is the building number of the former VDNKh exhibition hall on Yanki Kupaly st.
# This building is located at Yanki Kupaly, 27.
x = 27

# Step 2: Identify the number Y.
# Y is the number of a prominent former tram route on Kastryčnickaja st.
# Tram route number 4 used to pass through this street.
y = 4

# Step 3: Calculate X to the power of Y.
result = math.pow(x, y)

# Step 4: Print the final equation as requested.
# We cast the result to an integer because the result of a power operation on integers is an integer.
print(f"{x} ^ {y} = {int(result)}")