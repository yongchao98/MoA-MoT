def get_upper_group_order(v):
    """Returns the order of the upper ramification group G^v."""
    if v <= 1:
        return 8  # D4
    elif v <= 3:
        return 4  # C4
    else:
        return 1  # trivial group {e}

def phi(v):
    """Herbrand's phi-function to map upper to lower numbering."""
    if v <= 1:
        # For v in [0,1], d(phi)/dv = [D4:D4] = 1. phi(v) = v.
        return v
    elif v <= 3:
        # For v in (1,3], d(phi)/dv = [D4:C4] = 2.
        # phi(v) = phi(1) + integral from 1 to v of 2 du = 1 + 2*(v-1) = 2v-1.
        return 2 * v - 1
    else:
        # For v > 3, d(phi)/dv = [D4:{e}] = 8.
        # phi(v) = phi(3) + integral from 3 to v of 8 du = 5 + 8*(v-3) = 8v-19.
        return 8 * v - 19

def psi(i):
    """Inverse of phi-function, to map lower to upper numbering."""
    if i <= 1:
        # i = v
        return i
    elif i <= 5:
        # i = 2v-1 => v = (i+1)/2
        return (i + 1) / 2
    else:
        # i = 8v-19 => v = (i+19)/8
        return (i + 19) / 8

def get_lower_group_order(i):
    """Returns the order of the lower ramification group G_i."""
    v = psi(i)
    return get_upper_group_order(v)

def solve():
    """
    Finds the smallest integer t for which the lower ramification filtration is trivial.
    """
    t = 1
    while True:
        order = get_lower_group_order(t)
        # print(f"G_{t} has order {order}")
        if order == 1:
            print(f"The smallest integer t for which the lower filtration is trivial is: {t}")
            break
        t += 1

solve()