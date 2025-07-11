def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the ex-dividend share price in year 1 under the new policy.

    Args:
        q (float): The initial number of outstanding shares.
        E (float): The total market value of the firm's equity.
        d (float): The total dividend to be paid in year 1 under the original policy.
        g (float): The annual growth rate of the dividends.
    """
    # These are example values. You can replace them with your specific inputs.
    # q = 1000000.0
    # E = 50000000.0
    # d = 2000000.0
    # g = 0.05

    if q <= 0:
        print("Number of shares (q) must be positive.")
        return

    # Calculate the ex-dividend price p1 using the derived formula
    p1 = (E * (1 + g) - d) / q

    # Print the final equation with the variables substituted
    print("The formula for the ex-dividend price p1 is: (E * (1 + g) - d) / q")
    print("\nSubstituting the given values:")
    print(f"p1 = ({E} * (1 + {g}) - {d}) / {q}")
    print(f"p1 = ({E * (1 + g)} - {d}) / {q}")
    print(f"p1 = {E * (1 + g) - d} / {q}")
    print("\nThe calculated ex-dividend price per share in year 1 is:")
    print(f"p1 = {p1}")

# --- User Input Section ---
# Please define the values for your specific case here.
# For example:
# q_val = 1000000.0
# E_val = 50000000.0
# d_val = 2000000.0
# g_val = 0.05
# calculate_ex_dividend_price(q_val, E_val, d_val, g_val)

# To make the code runnable with example values, we can call the function with sample data.
print("Running calculation with example data:")
calculate_ex_dividend_price(q=1000000.0, E=50000000.0, d=2000000.0, g=0.05)