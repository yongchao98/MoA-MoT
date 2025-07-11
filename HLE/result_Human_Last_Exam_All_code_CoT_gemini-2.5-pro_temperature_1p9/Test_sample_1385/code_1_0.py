def solve():
    """
    Determines and prints the correct order of Old Russian enclitics.
    """
    # The list of enclitics in their correct hierarchical order
    ordered_enclitics = ["же", "бо", "бы", "мя", "еси"]

    # Print the ordered list as a comma-separated string
    result = ", ".join(ordered_enclitics)
    print(result)

    # Final answer in the specified format
    print(f"<<<{result}>>>")

solve()