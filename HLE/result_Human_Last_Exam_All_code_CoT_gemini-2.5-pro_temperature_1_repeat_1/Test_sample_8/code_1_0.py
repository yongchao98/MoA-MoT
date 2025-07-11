import sympy

# Let's define the symbolic variable for alpha
alpha = sympy.Symbol('alpha')

# 1. Initial conormal order
initial_order_val = 2
initial_order = initial_order_val + alpha
print(f"The initial function f is in the conormal space A^s(X), where the conormal order s is: {initial_order}")
print("-" * 30)

# 2. Order of the operators
# The wave operator Box_g is a second-order differential operator.
operator_order_m = 2
print(f"The differential operator P = Box_g - sigma^2 has order m = {operator_order_m}.")

# The resolvent R(sigma) = P^-1 is a pseudodifferential operator.
# Its order 'k' is the negative of the order of P.
resolvent_order_k = -operator_order_m
print(f"The resolvent R(sigma) is a pseudodifferential operator of order k = -m = {resolvent_order_k}.")
print("-" * 30)

# 3. Mapping property of pseudodifferential operators
print("An operator of order k applied to a function in A^s(X) results in a function in A^(s-k)(X).")
print("We need to calculate the new conormal order s' = s - k.")
print("-" * 30)

# 4. Calculation of the final conormal order
print("The final equation for the new order s' is:")
# The instruction is to output each number in the final equation.
# s' = (2 + alpha) - (-2)
final_order = initial_order - resolvent_order_k
print(f"s' = ({initial_order_val} + alpha) - ({resolvent_order_k})")
print(f"s' = {initial_order_val} + alpha + {abs(resolvent_order_k)}")
print(f"s' = {final_order}")
print("-" * 30)

# 5. Conclusion
print("Conclusion:")
print(f"The function R(sigma)f belongs to the conormal space A^s'(X) with s' = {final_order}.")
print(f"So, R(sigma)f is in A^({final_order})(X).")
