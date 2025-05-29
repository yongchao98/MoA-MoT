# Given values
HP = 5
HQ = 2

# Using the properties of the orthocenter and altitudes
# BP = HP + PC and AQ = HQ + QC
# We need to find (BP)(PC) - (AQ)(QC)

# Assume PC = x and QC = y
# Then BP = HP + x = 5 + x
# AQ = HQ + y = 2 + y

# We need to find (5 + x)x - (2 + y)y
# Since HP * HA = HQ * HB, we can assume HA = 2 and HB = 5 for simplicity
# This gives us x = 2 and y = 5

# Calculate the expression
x = 2
y = 5
BP = HP + x
AQ = HQ + y

result = (BP * x) - (AQ * y)
print(result)