import struct

def simulate_c_code_fix():
    """
    This script simulates the memory operations of the given C code
    to explain the instability and the fix.
    """
    print("--- Simulating the C Code ---")

    # In C: unsigned long long n = 1 << 10;
    # In Python, we simulate this with an 8-byte bytearray.
    # '<Q' means little-endian unsigned long long (8 bytes).
    n = 1024
    n_bytes = bytearray(struct.pack('<Q', n))
    print(f"Initial bytes of 'n' (value {n}): {list(n_bytes)}")

    # In C: process(&n); which does *p = 37; p[1] = 100;
    # We modify the first two bytes of our bytearray.
    n_bytes[0] = 37  # ASCII for '%'
    n_bytes[1] = 100 # ASCII for 'd'
    print(f"Bytes of 'n' after 'process()': {list(n_bytes)}")

    # The bytes now represent the null-terminated string "%d".
    # We find the first null byte to extract the string.
    try:
        end_of_string = n_bytes.index(0)
        format_string = n_bytes[:end_of_string].decode('ascii')
    except ValueError:
        format_string = n_bytes.decode('ascii')
    
    print(f"Interpreted as a string, this is: '{format_string}'")
    print("\n--- Analyzing printf() Behavior ---")

    # Original unstable code: printf(s, s);
    # This prints the memory address of 'n', which changes on each run due to ASLR.
    print("Original code prints the memory address of 'n'. This output is UNSTABLE.")
    
    # Fixed stable code: printf(s, *s);
    # *s dereferences the pointer, getting the value of the first byte.
    value_to_print = n_bytes[0]
    print("\nFixed code prints the value pointed to by 's' (the first byte).")
    print(f"Format String: '{format_string}'")
    print(f"Value to Print: {value_to_print}")
    
    # Python's equivalent of C's printf("%d", 37)
    final_output = format_string % value_to_print
    
    print("\n--- Final Stable Output ---")
    print(final_output)

simulate_c_code_fix()