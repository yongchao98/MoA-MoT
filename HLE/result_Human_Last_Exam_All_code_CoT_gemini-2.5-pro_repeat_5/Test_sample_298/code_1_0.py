def format_group(name, order, count):
    """Formats a component of the group string."""
    if count == 0:
        return ""
    if name == "Z":
        if count == 1:
            return "Z"
        else:
            return "Z^" + str(count)
    if order == 0: # Should not happen
        return ""
    group_str = f"{name}/{order}{name}"
    if count > 1:
        return "+".join([group_str] * count)
    else:
        return group_str

def get_cohomology_groups(k):
    """
    Calculates the list of cohomology groups H^i(M(k))
    based on known results for the homology of the braid group of the projective plane.
    This function implements the logic derived in the text.
    """
    if k == 0:
        return "[Z]"
    if k == 1:
        # M(1) is homotopy equivalent to S^1
        return "[Z, Z]"

    # Homology groups H_q(B_k(RP^2)) for k>=4
    # H_q = 0 for q >= k
    H = {0: ("Z", 1)}
    if k >= 2:
        H[1] = [("Z", 1), ("Z/2Z", 1)]
    if k >= 4:
        H[2] = [("Z/2Z", 1)]
        for j in range(2, (k // 2) + 2):
            if 2 * j - 1 < k:
                H[2 * j - 1] = [("Z/2Z", j)]
            if 2 * j < k:
                H[2 * j] = [("Z/2Z", j)]

    # Compute Cohomology H^q using Universal Coefficient Theorem
    # H^q = Hom(H_q, Z) + Ext(H_{q-1}, Z)
    H_coho = {}

    # H^0
    H_coho[0] = [("Z", 1)]
    
    # H^1
    H_coho[1] = [("Z", H[1][0][1])] # Hom(Z,Z) part from H_1

    # H^q for q >= 2
    for q in range(2, k + 1):
        ext_parts = []
        if q - 1 in H:
            h_prev = H[q-1]
            for group, count in h_prev:
                if "Z/2Z" in group:
                    ext_parts.append(("Z/2Z", count))
        if ext_parts:
            H_coho[q] = ext_parts

    # Format the final output string
    result = []
    for i in range(k + 1):
        if i in H_coho:
            parts = []
            for group, count in H_coho[i]:
                if group == "Z":
                    parts.append(format_group("Z", 0, count))
                elif "Z/2Z" in group:
                    parts.append(format_group("Z", 2, count))
            if not parts:
                result.append("0")
            else:
                result.append("+".join(p for p in parts if p))
        else:
            result.append("0")

    # Remove trailing zeros
    last_nonzero = -1
    for i in range(len(result) - 1, -1, -1):
        if result[i] != "0":
            last_nonzero = i
            break
    
    return "[" + ", ".join(result[:last_nonzero+1]) + "]"

k = 7
cohomology_list_str = get_cohomology_groups(k)
print(cohomology_list_str)