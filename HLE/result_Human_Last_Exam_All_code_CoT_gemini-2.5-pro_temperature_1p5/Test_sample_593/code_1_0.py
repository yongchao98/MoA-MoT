def calculate_treewidth_bound(t_H, t_G, k):
    """
    Calculates the tight upper bound on the treewidth of the union of two graphs.

    Args:
        t_H (int): The treewidth of graph H.
        t_G (int): The treewidth of graph G.
        k (int): The number of shared vertices between H and G.

    Returns:
        int: The tight upper bound for the treewidth of the combined graph F.
    """
    # The tight upper bound is given by max(t_H, t_G) + k.
    upper_bound = max(t_H, t_G) + k
    return upper_bound

def main():
    """
    Main function to demonstrate the calculation with example values.
    """
    # Example values. You can change these to test other scenarios.
    t_H = 7
    t_G = 5
    k = 4

    print("This script calculates the tight upper bound for the treewidth of a graph F,")
    print("which is the union of two graphs H and G sharing k vertices.")
    print("The formula for the bound is: tw(F) <= max(t_H, t_G) + k\n")

    # Calculate the bound
    bound = calculate_treewidth_bound(t_H, t_G, k)

    # Print the equation with the specific numbers
    print(f"Given values:")
    print(f"  Treewidth of H (t_H): {t_H}")
    print(f"  Treewidth of G (t_G): {t_G}")
    print(f"  Number of shared vertices (k): {k}\n")
    print(f"The equation for the bound is:")
    print(f"tw(F) <= max({t_H}, {t_G}) + {k}")
    print(f"tw(F) <= {max(t_H, t_G)} + {k}")
    print(f"tw(F) <= {bound}\n")
    print(f"The tight upper bound is: {bound}")


if __name__ == "__main__":
    main()
