def generate_ita_matrix_code():
    """
    This function generates the ITA Matrix extension code for finding
    direct Delta flights eligible for Platinum Companion Certificates,
    following all specified constraints.
    """

    # 1. Define the eligible fare classes for Platinum Companion Certificates.
    # These must be uppercase.
    fare_classes = ['L', 'U', 'T', 'X', 'V']

    # 2. To get the lexicographically highest string of minimal length,
    # sort the fare classes in descending (reverse) alphabetical order.
    fare_classes.sort(reverse=True)

    # 3. Join the sorted classes with the '|' (OR) separator.
    class_string = '|'.join(fare_classes)

    # 4. Construct the final extension code.
    # 'f' is the command for fare rules.
    # 'c' is the shortest alias for booking-class.
    # Using lowercase and no spaces minimizes length and adheres to instructions.
    extension_code = f"f c={class_string}"

    print(extension_code)

generate_ita_matrix_code()