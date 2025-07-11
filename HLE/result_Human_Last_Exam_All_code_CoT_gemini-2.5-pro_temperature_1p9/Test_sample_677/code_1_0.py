def find_factors(target_area):
    """
    Finds if a target_area can be formed by multiplying two integers n, m >= 4.
    If it can, prints the equation and returns True. Otherwise returns False.
    """
    found = False
    # We only need to check for n up to sqrt(target_area)
    for n in range(4, int(target_area**0.5) + 1):
        if target_area % n == 0:
            m = target_area // n
            if m >= 4:
                print(f"The area {target_area} can be formed by the rectangle: {n} x {m}")
                found = True
                # We can break if we just need one example, but let's find all.
    return found

def main():
    """
    Main function to check which board sizes lead to a planar graph.
    """
    print("Checking for largest area of a planar graph...")
    print("Known non-planar graphs have areas nm >= 35.")
    print("We search for the largest possible area strictly less than 35.")

    supremum_candidate = 34
    print(f"\nChecking if area {supremum_candidate} is possible for n, m >= 4:")
    if not find_factors(supremum_candidate):
        print(f"No integer rectangle with n, m >= 4 has an area of {supremum_candidate}.")
        # Let's check the next candidate down
        supremum_candidate = 33
        print(f"\nChecking if area {supremum_candidate} is possible for n, m >= 4:")
        if not find_factors(supremum_candidate):
            print(f"No integer rectangle with n, m >= 4 has an area of {supremum_candidate}.")
            # And the next
            supremum_candidate = 32
            print(f"\nChecking if area {supremum_candidate} is possible for n, m >= 4:")
            if find_factors(supremum_candidate):
                print(f"The graph for the {supremum_candidate} area board is planar.")
                print(f"\nThe largest area nm for which the graph is planar, given that areas >= 35 are not, is {supremum_candidate}.")

if __name__ == "__main__":
    main()
