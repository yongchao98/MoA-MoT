# The user wants to multiply two 100-digit numbers, p and q.
# Since specific values for p and q were not provided, we will use
# two example 100-digit numbers for the demonstration.
# Python has native support for arbitrarily large integers, which simplifies the multiplication.

# Input numbers p and q
p_str = "1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890"
q_str = "9876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210"

# Convert the string representations to integer types for calculation
p = int(p_str)
q = int(q_str)

# Perform the multiplication
o = p * q

# Output the final equation as requested.
# The numbers are very large, so we print them on separate lines for clarity.

print("p:")
print(p)
print("\nq:")
print(q)
print("\no = p * q:")
print(o)

# The final part of the task is to calculate the minimized total memory use (m)
# for p, q, and o in 'D' (decimal digits) for the Wuxing architecture.
#
# 1. Memory for p: p has 100 digits. This requires 100 D of storage.
# 2. Memory for q: q has 100 digits. This requires 100 D of storage.
# 3. Memory for o: The product of two 100-digit numbers has a maximum of
#    100 + 100 = 200 digits. This requires 200 D of storage.
#
# Total minimized memory m = 100 D (for p) + 100 D (for q) + 200 D (for o) = 400 D.