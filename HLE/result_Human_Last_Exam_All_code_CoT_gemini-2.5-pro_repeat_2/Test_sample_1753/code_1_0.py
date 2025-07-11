import math

# The derived relationship between the arc length L and the constant a is L = 3 * a^(2/3).
# We are given that L = 3/2.

L = 3/2
print(f"The relationship between arc length L and the constant a is: L = 3 * a^(2/3)")
print(f"Given L = 3/2, we set up the equation:")
# The numbers in this equation are 3, 2, 3, 3, 2
print(f"3 * a^(2/3) = 3/2")

# To solve for a, first we solve for a^(2/3)
# a^(2/3) = (3/2) / 3 = 1/2
a_pow_2_3 = (3/2) / 3
print(f"\nStep 1: Solve for a^(2/3)")
print(f"a^(2/3) = (3/2) / 3")
# The numbers in this equation are 2, 3, 1, 2
print(f"a^(2/3) = 1/2")

# Now, we solve for a by raising both sides to the power of 3/2
# a = (1/2)^(3/2)
a = (1/2)**(3/2)
print(f"\nStep 2: Solve for a by raising both sides to the power of 3/2")
print(f"a = (1/2)^(3/2)")

print(f"\nThe final value for a is {a}")
print(f"This is equivalent to 1/(2*sqrt(2)) or sqrt(2)/4.")

# As requested, here are the numbers in the final equation a = (1/2)^(3/2)
print("\nThe numbers in the final equation are:")
print(f"Base numerator: 1")
print(f"Base denominator: 2")
print(f"Exponent numerator: 3")
print(f"Exponent denominator: 2")