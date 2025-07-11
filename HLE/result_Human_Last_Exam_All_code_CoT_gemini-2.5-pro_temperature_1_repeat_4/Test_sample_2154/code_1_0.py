import math

def u_r(n):
  """
  Calculates the minimal order of the Picard-Fuchs differential equation
  for the given Hamiltonian, based on the formula from OEIS A115252.
  """
  return math.floor((n + 1) / 3)

# We need to find the set of values for n from 3 to 12.
results_list = []
for n in range(3, 13):
  value = u_r(n)
  results_list.append(value)

# The problem asks to "Find {u_r(3), u_r(4), ..., u_r(12)}".
# The final result is the set of these numbers.
# The following print statements fulfill the requirement to "output each number".
print("The values for the set {u_r(3), u_r(4), ..., u_r(12)} are:")
print("{", end="")
for i, result in enumerate(results_list):
    if i > 0:
        print(", ", end="")
    print(result, end="")
print("}")

# For clarity, here is each individual calculation:
print("\nDetailed calculation:")
for n in range(3, 13):
    print(f"u_r({n}) = floor(({n}+1)/3) = {(n+1)//3}")
