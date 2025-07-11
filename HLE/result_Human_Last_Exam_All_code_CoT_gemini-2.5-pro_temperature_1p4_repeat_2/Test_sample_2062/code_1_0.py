import math

def count_similarity_classes(n, q):
  """
  Calculates the number of similarity classes for diagonal n x n matrices
  over a field of size q. This is equivalent to C(n + q - 1, n).
  """
  if n < 0 or q < 1:
    return 0
  return math.comb(n + q - 1, n)

def demonstrate_growth():
  """
  Demonstrates that the number of similarity classes grows polynomially,
  not exponentially, by comparing it with 2^n.
  """
  q_fixed = 4  # Example: a field with 4 elements
  base_exponential = 2
  
  print(f"Demonstrating growth for fixed q = {q_fixed}:")
  print("-" * 50)
  print(f"{'n':>5} | {'Num. Classes (Polynomial)':>25} | {'Exponential (' + str(base_exponential) + '^n)':>20}")
  print("-" * 50)
  
  for n in range(1, 21, 2):
    poly_growth = count_similarity_classes(n, q_fixed)
    exp_growth = base_exponential**n
    print(f"{n:>5} | {poly_growth:>25,} | {exp_growth:>20,}")
  print("-" * 50)
  print("As you can see, the number of classes grows much slower than the exponential function.\n")

# Run the demonstration
demonstrate_growth()

# The final answers to the questions
answer_a = "Yes"
answer_b = 1
answer_c = "No"

# Format and print the final answer as requested
# For (b), we output the number 1 from the variable 'answer_b'.
# The final "equation" is the formatted answer string.
final_answer_string = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
print("Final Answer:")
print(final_answer_string)
print("<<<" + f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}" + ">>>")