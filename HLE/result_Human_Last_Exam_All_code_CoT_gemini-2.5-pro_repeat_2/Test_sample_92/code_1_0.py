def solve():
    """
    Calculates the probability based on the linear potential solution to the gambler's ruin problem.
    """
    start_bin = 0
    melt_bin = 2024
    escape_bin = 2025

    # The probability h(n) of reaching escape_bin before melt_bin, starting from n,
    # for this type of symmetric random walk is given by the formula:
    # h(n) = (n - melt_bin) / (escape_bin - melt_bin)
    
    # We need to calculate h(start_bin)
    numerator = start_bin - melt_bin
    denominator = escape_bin - melt_bin
    
    result = numerator / denominator

    # The problem asks for the final equation, so we print the components.
    print(f"The probability is calculated using the formula:")
    print(f"P = (start_bin - melt_bin) / (escape_bin - melt_bin)")
    print(f"P = ({start_bin} - {melt_bin}) / ({escape_bin} - {melt_bin})")
    print(f"P = {numerator} / {denominator}")
    print(f"P = {result}")

solve()