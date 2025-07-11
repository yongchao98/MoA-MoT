def F(x):
  """
  Represents the predicate 'is greater than 10'.
  This is our concept 'F'.
  """
  return x > 10

# 1. Understanding the proposition Fa
# Here, 'a' is an individual concept, represented by the number 20.
a = 20

# We can understand and evaluate 'Fa' by combining the concepts.
understanding_Fa = F(a)
print("Understanding the specific proposition 'Fa':")
print(f"The proposition is F(a), where F is 'x > 10' and a = {a}.")
print(f"The equation is: {a} > 10 = {understanding_Fa}")
print("-" * 20)

# 2. Understanding the proposition ∀x Fx (For all x, Fx is true)
# This requires understanding F (which we have) and universal quantification (∀).
# In Python, `all()` is a perfect model for universal quantification.
# Let's define a domain of discourse for 'x'.
domain_of_x = [15, 30, 55, 100]

# If we have the concept 'F' and understand quantification ('all'),
# we can construct and evaluate the universal proposition.
understanding_forall_x_Fx = all(F(x) for x in domain_of_x)

print("Understanding the universal proposition '∀x Fx':")
print(f"The proposition is ∀x Fx, over the domain {domain_of_x}.")
# We will show each number being evaluated in the final 'equation'.
evaluations = [f"{x} > 10 = {F(x)}" for x in domain_of_x]
print(f"The 'equation' involves checking each number: all([{', '.join(evaluations)}])")
print(f"The final result is: {understanding_forall_x_Fx}")
print("-" * 20)

# Example of a failing case for clarity
domain_with_failure = [15, 30, 5, 55]
understanding_fail_case = all(F(x) for x in domain_with_failure)
print("A failing case for '∀x Fx':")
print(f"Domain: {domain_with_failure}")
evaluations_fail = [f"{x} > 10 = {F(x)}" for x in domain_with_failure]
print(f"The 'equation': all([{', '.join(evaluations_fail)}])")
print(f"The final result is: {understanding_fail_case}")
