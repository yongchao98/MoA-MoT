# Plan:
# 1. The initial function f belongs to the conormal space A^{2+alpha}(X). The regularity index is s = 2 + alpha.
# 2. The resolvent R(sigma) is associated with the wave operator, which is a second-order differential operator. The resolvent is a pseudodifferential operator of order -2.
# 3. Applying an elliptic pseudodifferential operator of order m to a function in a conormal space A^s results in a function in the space A^{s-m}.
# 4. We have s = 2 + alpha and m = -2. The new regularity index is s' = s - m = (2 + alpha) - (-2) = 4 + alpha.
# 5. The code will calculate this new index and print the name of the resulting space.

# The variable alpha is a real number in (0, 1). We'll represent it as a string.
alpha_str = "alpha"

# The integer part of the initial regularity index s = 2 + alpha.
initial_order_part = 2

# The resolvent R(sigma) is an operator of order -2. This is the key property.
operator_order = -2

# The new regularity index is s_new = s - operator_order.
# The integer part of the new regularity index is calculated here.
final_order_part = initial_order_part - operator_order

# The numbers in the final equation are the new order part (4) and alpha.
# The final equation describes the space A^{4+alpha}(X).
# The code below prints this result.
print(f"The resolvent R(σ)f belongs to the conormal space A^({final_order_part}+{alpha_str})(X).")
