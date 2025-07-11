def generate_sequence_s4(num_terms):
    """
    Generates and prints the terms of sequence S4 based on the deduced rule.
    The rule is s[n] = s[n-2] + s[n-3].
    The base cases are s[1]=1, s[2]=1, s[3]=2 from the provided sequence.
    """
    s = {1: 1, 2: 1, 3: 2}

    print("s[1] = 1")
    print("s[2] = 1")
    print("s[3] = 2")

    for n in range(4, num_terms + 1):
        s_n_minus_2 = s[n - 2]
        s_n_minus_3 = s[n - 3]
        s[n] = s_n_minus_2 + s_n_minus_3
        print(f"s[{n}] = s[{n - 2}] + s[{n - 3}] = {s_n_minus_2} + {s_n_minus_3} = {s[n]}")

# The provided sequence S4 has 42 terms.
# We will generate the same number of terms using the deduced rule.
generate_sequence_s4(42)
