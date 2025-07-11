import math

def get_homotopy_group_size_str(k):
    """Returns the size of the k-th homotopy group of SO(3) as a string for display."""
    if k < 0:
        return "1"  # By convention, pi_k for k<0 is trivial.
    sizes = {
        0: "1",       # |pi_0(SO(3))| = 1 (connected)
        1: "2",       # |pi_1(SO(3))| = |Z_2| = 2
        2: "1",       # |pi_2(SO(3))| = |{0}| = 1
        3: "|Z|",     # |pi_3(SO(3))| = |Z| = infinity
        4: "2",       # |pi_4(SO(3))| = |Z_2| = 2
    }
    return sizes[k]

def get_homotopy_group_size_val(k):
    """Returns the numerical size of the k-th homotopy group of SO(3)."""
    if k < 0:
        return 1
    sizes = {
        0: 1,
        1: 2,
        2: 1,
        3: math.inf,
        4: 2,
    }
    return sizes[k]

def get_group_name(k):
    """Returns the string name of the k-th homotopy group."""
    return f"pi_{k}(SO(3))"

print("The number of circulons for G=SO(3) in d+1 dimensions is calculated as follows:")
print("-" * 40)

# Case d=1
d = 1
print(f"For d={d}:")
print("A circle cannot be embedded in R^1. This is a degenerate case.")
print("By convention, there is only 1 trivial configuration.")
print("Number of circulons = 1")
print("-" * 40)

# Case d=2
d = 2
print(f"For d={d}:")
print("The classification is given by the homotopy classes of maps [R^2 \\ S^1, SO(3)],")
print("which corresponds to |pi_0(SO(3))| * |pi_1(SO(3))|.")
p0_size_str = get_homotopy_group_size_str(0)
p1_size_str = get_homotopy_group_size_str(1)
p0_size_val = get_homotopy_group_size_val(0)
p1_size_val = get_homotopy_group_size_val(1)
result = p0_size_val * p1_size_val
print(f"Number of circulons = |{get_group_name(0)}| * |{get_group_name(1)}| = {p0_size_str} * {p1_size_str} = {result}")
print("-" * 40)

# Cases d=3 to d=6
for d in range(3, 7):
    k = d - 2
    print(f"For d={d}:")
    print(f"The classification is given by the homotopy classes of maps [R^{d} \\ S^1, SO(3)],")
    print(f"which corresponds to |pi_{k}(SO(3))| * |pi_1(SO(3))|.")
    pk_size_str = get_homotopy_group_size_str(k)
    p1_size_str = get_homotopy_group_size_str(1)
    pk_size_val = get_homotopy_group_size_val(k)
    p1_size_val = get_homotopy_group_size_val(1)
    
    result = pk_size_val * p1_size_val
    if math.isinf(result):
        result_str = "infinity"
    else:
        result_str = str(int(result))

    print(f"Number of circulons = |{get_group_name(k)}| * |{get_group_name(1)}| = {pk_size_str} * {p1_size_str} = {result_str}")
    print("-" * 40)