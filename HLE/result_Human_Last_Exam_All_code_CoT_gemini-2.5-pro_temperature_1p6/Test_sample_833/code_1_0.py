import sys

# We analyze the expression at the shock front.
# Let's calculate the value on the right side of the shock (x -> 0+).
u_right = 0
# For the shock u(x) = H(-x), u_x = -delta(x).
# u_bar = K * u_x = -1/2 * exp(-|x|). At x=0, u_bar = -1/2.
u_bar_val = -1/2
# u_bar_x = 1/2 * sgn(x) * exp(-|x|). At x -> 0+, u_bar_x -> 1/2.
u_bar_x_val_right = 1/2

# The expression is (1-2u)*u_bar*u_bar_x
a_right = (1 - 2 * u_right) * u_bar_val * u_bar_x_val_right
print("Calculation for the right side of the shock (u=0):")
print(f"a = (1 - 2 * {u_right}) * ({u_bar_val}) * ({u_bar_x_val_right})")
print(f"a = {a_right}")
print("-" * 20)


# Let's calculate the value on the left side of the shock (x -> 0-).
u_left = 1
# u_bar is continuous, so u_bar = -1/2.
# u_bar_x at x -> 0- is -1/2.
u_bar_x_val_left = -1/2

# The expression is (1-2u)*u_bar*u_bar_x
a_left = (1 - 2 * u_left) * u_bar_val * u_bar_x_val_left
print("Calculation for the left side of the shock (u=1):")
print(f"a = (1 - 2 * {u_left}) * ({u_bar_val}) * ({u_bar_x_val_left})")
print(f"a = {a_left}")
print("-" * 20)

print(f"The lower bound is {a_right}.")

# Final Answer
sys.stdout.write("<<<-0.25>>>")