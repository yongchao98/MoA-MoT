import sys

# In a real script, n might be a user input, e.g., n = int(sys.argv[1])
# For demonstration, we'll use a fixed value.
n = 4

print(f"We are considering n = {n}.")
print("The set A is composed of n-tuples (x_1, ..., x_n) where each x_i is a rational cube.")
print("The condition for a tuple to be in A is the conjunction of n independent conditions:")
print("x_1 is a cube, AND x_2 is a cube, ..., AND x_n is a cube.")

print("\nEach condition 'x_i is a cube' can be expressed using one auxiliary variable y_i:")
print("∃ y_i in Q such that x_i - y_i^3 = 0.")

print("\nTo combine these n conditions into a single polynomial equation, we sum their squares.")
print("This works because for rational numbers, a sum of squares is zero if and only if each term is zero.")
print("This construction requires n auxiliary variables y_1, ..., y_n.")
print(f"This shows that m <= {n}.")

print("\nBecause the n conditions are independent, we need at least n auxiliary variables to act as 'witnesses' for each condition.")
print(f"This implies that m >= {n}.")
print(f"\nThus, the smallest possible value for m is exactly {n}.")
print("-" * 20)

m = n
print(f"The final answer is m = {m}")

print("\nThe polynomial equation F=0 that defines the set A is built as follows:")
equation_parts = []
for i in range(1, n + 1):
    equation_parts.append(f"(X{i} - Y{i}^{3})^{2}")

final_equation = " + ".join(equation_parts) + " = 0"
print(final_equation)

print("\nTo comply with the prompt, we list the numbers present in this equation:")
# The numbers are coefficients and exponents.
# They are independent of n.
print(" - The number 1 (as the implicit coefficient of X_i and Y_i^3)")
print(" - The number 3 (as the exponent of each Y_i variable)")
print(" - The number 2 (as the exponent of each squared term)")
print(" - The number 0 (as the right-hand side of the equation)")
