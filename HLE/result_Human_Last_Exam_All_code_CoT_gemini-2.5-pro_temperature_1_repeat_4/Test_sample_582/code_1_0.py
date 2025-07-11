import sys

def calculate_ex_dividend_price(q, E, d, g):
    """
    Calculates the per-share ex-dividend price in year 1 under the new policy.

    Args:
        q (float): The current number of outstanding shares.
        E (float): The total market value of the firm's equity.
        d (float): The total dividends to be distributed in year 1 under the old policy.
        g (float): The annual growth rate of the dividends.
    """
    # Check for division by zero
    if q == 0:
        print("Error: The number of shares (q) cannot be zero.")
        return

    # Numerator of the formula: E * (1 + g) - d
    numerator = E * (1 + g) - d
    # Denominator of the formula: q
    denominator = q

    # Calculate the per-share ex-dividend price in year 1
    p1 = numerator / denominator

    # Print the calculation steps
    print("The formula for the per-share ex-dividend price (p1) in year 1 is:")
    print("p1 = (E * (1 + g) - d) / q")
    print("\nSubstituting the given values:")
    print(f"p1 = ({E} * (1 + {g}) - {d}) / {q}")
    print(f"p1 = ({E * (1 + g)} - {d}) / {q}")
    print(f"p1 = {numerator} / {q}")
    print(f"p1 = {p1}")

    # The final answer in the requested format will be the symbolic expression,
    # but the script prints the numerical result for the user.

if __name__ == '__main__':
    # --- User-defined variables ---
    # You can change these values to match your specific case.
    # q: current number of outstanding shares
    q = 1000000.0
    # E: total market value of equity
    E = 50000000.0
    # d: total dividends in year 1 (original policy)
    d = 2000000.0
    # g: annual dividend growth rate (e.g., 0.05 for 5%)
    g = 0.05

    # Check if arguments are provided via command line
    if len(sys.argv) == 5:
        try:
            q = float(sys.argv[1])
            E = float(sys.argv[2])
            d = float(sys.argv[3])
            g = float(sys.argv[4])
        except ValueError:
            print("Invalid arguments. Please provide four numbers for q, E, d, and g.")
            print("Usage: python your_script_name.py <q> <E> <d> <g>")
            sys.exit(1)
    
    # Execute the calculation
    calculate_ex_dividend_price(q, E, d, g)
