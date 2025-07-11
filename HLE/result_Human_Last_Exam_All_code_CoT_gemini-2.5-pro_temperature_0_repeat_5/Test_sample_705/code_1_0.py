import itertools
import sys

def get_C_set_size(n, num_vals_per_component):
    """
    Calculates the size of the set C(sigma), where sigma contains
    num_vals_per_component values for each of the n components.
    This represents the number of states that must be simulated in one
    step of the relaxed simulation.
    """
    size = 1
    for _ in range(n):
        size *= num_vals_per_component
    return size

def get_memory_usage_bytes(num_states, n):
    """
    Estimates the memory to store a list of states, assuming each state
    component takes 8 bytes (e.g., a pointer or a 64-bit number).
    """
    # Memory for the list of pointers to state tuples
    list_overhead = sys.getsizeof([])
    pointer_size = 8 # 8 bytes per pointer in 64-bit Python
    list_pointers_mem = num_states * pointer_size

    # Memory for the state tuples themselves
    tuple_overhead = sys.getsizeof(tuple())
    tuple_elements_mem = n * 8 # 8 bytes per element
    states_mem = num_states * (tuple_overhead + tuple_elements_mem)

    return list_overhead + list_pointers_mem + states_mem

def format_bytes(byte_count):
    """Formats a byte count into KB, MB, GB, etc."""
    if byte_count is None:
        return "N/A"
    power = 1024
    n = 0
    power_labels = {0: '', 1: 'KB', 2: 'MB', 3: 'GB', 4: 'TB'}
    while byte_count >= power and n < len(power_labels):
        byte_count /= power
        n += 1
    return f"{byte_count:.2f} {power_labels[n]}"


print("Demonstrating the potential exponential growth of memory space for computation.")
print("The intermediate set C(sigma) often needs to be generated and stored in memory.")
print("Let's assume each component Vk can take one of 2 values, and our current")
print("set of values sigma contains both values for each of the n components.")
print("-" * 70)
print(f"{'n (components)':<18} | {'|C(sigma)| (num states)':<25} | {'Estimated Memory Usage':<25}")
print("-" * 70)

for n in range(1, 31):
    # Calculate the number of states in the intermediate set C(sigma)
    num_states = get_C_set_size(n, 2)
    
    # Estimate memory usage for storing these states
    # Stop if it gets too large to avoid large number issues
    mem_usage = None
    if n < 30: # Prevent overflow for very large n
        mem_usage = get_memory_usage_bytes(num_states, n)

    print(f"{n:<18} | {num_states:<25} | {format_bytes(mem_usage):<25}")
    if n == 10:
        print("...")
    if n > 10 and n < 28:
        continue
