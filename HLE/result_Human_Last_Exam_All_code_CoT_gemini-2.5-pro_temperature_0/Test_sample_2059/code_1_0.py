# We want to compute the 20th term of a sequence u_n defined by:
# u_n = 6 * u_{n-1} - 2 * u_{n-2}
# with initial conditions u_0 = 1 and u_1 = 4.

# Initialize the first two terms of the sequence
u_prev = 1  # This is u_0
u_curr = 4  # This is u_1

# We need to calculate up to u_20. The loop will run from n=2 to n=19.
# After the loop, u_curr will hold u_19 and u_prev will hold u_18.
for n in range(2, 20):
    u_next = 6 * u_curr - 2 * u_prev
    u_prev = u_curr
    u_curr = u_next

# Now, u_curr is u_19 and u_prev is u_18
u_18 = u_prev
u_19 = u_curr

# The final calculation for u_20
u_20 = 6 * u_19 - 2 * u_18

# As requested, output each number in the final equation
print(f"The sum of squares is calculated from the recurrence u_n = 6 * u_{n-1} - 2 * u_{n-2}.")
print(f"The final step is u_20 = 6 * u_19 - 2 * u_18.")
print(f"Plugging in the values:")
print(f"6 * {u_19} - 2 * {u_18} = {u_20}")