from fractions import Fraction

# Step 1: Define the values for e(M_3) and e(H_3)
# e(M_3) = -B_6/6 = -(1/42)/6 = -1/252
e_M3 = Fraction(-1, 252)

# e(H_3) = -1/120
e_H3 = Fraction(-1, 120)

# Step 2: Calculate the orbifold Euler characteristic of the quotient stack [U/G]
# e([U/G]) = e(M_3) - e(H_3)
result = e_M3 - e_H3

# Step 3: Print the final equation with all the numbers
print("The orbifold Euler characteristic of the quotient stack [U/G] is calculated as follows:")
print(f"e([U/G]) = e(M_3) - e(H_3)")
print(f"e(M_3) = {e_M3}")
print(f"e(H_3) = {e_H3}")
print(f"e([U/G]) = {e_M3} - ({e_H3}) = {e_M3 + e_H3}")
print(f"e([U/G]) = {result}")

# Final Answer
final_answer_str = f"<<<{result}>>>"