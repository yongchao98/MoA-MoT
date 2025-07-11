# The order of the cyclic group C is a prime number p.
p = 31

# Based on the analysis, the group G is trivial.
# This means the set of extensions E has only one element, C_p itself.
# The sum is therefore the order of the outer automorphism group of C_p.
# o(C_p) = |Out(C_p)|
# For an abelian group, Out(H) = Aut(H).
# For a cyclic group C_p where p is prime, |Aut(C_p)| = phi(p) = p - 1.

# These are the numbers in the final equation:
# sum = o(C_31) = phi(31) = 31 - 1 = 30
val_p = 31
subtrahend = 1
result = val_p - subtrahend

# Print the numbers that form the final equation
print(f"The equation for the sum is: phi({val_p}) = {val_p} - {subtrahend} = {result}")

# Print the final result
print("The computed sum is:")
print(result)
