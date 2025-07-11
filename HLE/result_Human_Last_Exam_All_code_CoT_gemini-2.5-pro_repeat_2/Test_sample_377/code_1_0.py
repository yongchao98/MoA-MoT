import math

# The problem reduces to finding the number of conjugacy classes
# of the extraspecial group S = 3^{1+2}_+.

# An extraspecial group of order p^(2n+1) has a known formula for its
# number of conjugacy classes.
# For S = 3^(1+2), we have p=3 and 2n=2, so n=1.

p = 3
n = 1

# The number of conjugacy classes is given by the formula: p^(2n) + p - 1.
# We can also calculate this from first principles for n=1:
# The center Z(S) has p elements, forming p classes of size 1.
# The number of non-central elements is p^(2n+1) - p.
# For n=1, each non-central class has size p.
# So, the number of non-central classes is (p^3 - p) / p = p^2 - 1.
# Total classes = p + (p^2 - 1) = p^2 + p - 1.

# Let's use the formula to calculate the result.
num_classes = p**(2 * n) + p - 1

print("The number of blocks of kG is equal to the number of conjugacy classes of S.")
print(f"S is an extraspecial group of order p^(2n+1) with p={p}, n={n}.")
print(f"The number of conjugacy classes is given by the formula: p^(2n) + p - 1")
print(f"Substituting the values: {p}^(2*{n}) + {p} - 1")
print(f"= {p**(2*n)} + {p} - 1")
print(f"= {p**(2*n) + p} - 1")
print(f"= {num_classes}")

# The number of blocks does kG have is the final calculated number.
print(f"\nThus, the number of blocks kG has is {num_classes}.")
